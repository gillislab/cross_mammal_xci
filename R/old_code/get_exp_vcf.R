library(stringr)

BASES = c("A", "C", "G", "T")
POS_IDX = c(1)
FW_STRAND_IDX = 2:5
RV_STRAND_IDX = 9:12

VCF_POS_COL = 'V2'
VCF_REF_COL = 'V4'
VCF_ALT_COL = 'V5'

EXPRESSION_THRESHOLD = 10

get_exp_vcf = function(vcf_file_path, wig_file_path) {
  
  #Read in the vcf and wig files
  vcf = try(read.table(vcf_file_path))
  if(is.null(dim(vcf))){return(NA)}
  #filtering out Variants with , in them, they are multiple alleles missed previously
  vcf = vcf[!grepl(',', vcf[,VCF_REF_COL]),] #Ref alleles with commas
  vcf = vcf[!grepl(',', vcf[,VCF_ALT_COL]),] #Alt alleles with commas
  
  vcf = vcf[c(VCF_POS_COL, VCF_REF_COL, VCF_ALT_COL)]
  colnames(vcf) = c('Pos', 'ref', 'alt')
  
  wig = try(read.table(wig_file_path, skip = 3))
  if(is.null(dim(wig))){return(NA)}
  
  # Combine negative and positive strand counts
  wig_counts = wig[,FW_STRAND_IDX] + wig[,RV_STRAND_IDX]       
  wig_counts = cbind(wig[,POS_IDX], wig_counts)
  colnames(wig_counts) = c(c('Pos'), BASES)
  
  merged_df = merge(x=wig_counts, y=vcf, by="Pos")
  
  #Delineate the number of alleles at a position
  n_allele = rowSums(merged_df[,BASES]>0)
  # Check and keep positions with only two alleles
  biallelic_mask =  n_allele == 2
  merged_df = merged_df[biallelic_mask,]
  
  expression_mask = rowSums(merged_df[,BASES] >= EXPRESSION_THRESHOLD) == 2
  merged_df = merged_df[expression_mask,]
  
  merged_df['ref_counts'] = apply(merged_df, 1, function (x) {as.numeric(x[x['ref']])})
  merged_df['alt_counts'] = apply(merged_df, 1, function (x) {as.numeric(x[x['alt']])})
  
  merged_df =  merged_df[c('Pos','ref_counts','alt_counts','A','C','G','T', 'ref','alt')]
  
  #Any remaining SNP with fewer reads than EXPRESSION_THRESHOLD for the reference or alternate allele indicates a disagreement between the VCF 
  #identified alt and ref allele and the bases with read pileups
  #Just filter them out
  merged_df = merged_df[pmin(merged_df$ref_counts, merged_df$alt_counts) >= EXPRESSION_THRESHOLD,]
  return(merged_df)
}
