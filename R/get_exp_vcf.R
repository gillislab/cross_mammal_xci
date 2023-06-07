
#Function to stich together the vcf and wig output from the GATK pipeline to get the final variant and allelic-expression data

get_exp_vcf = function(vcf_file_path, wig_file_path){
  
  #Read in the vcf and wig files
  vcf = read.table(vcf_file_path)
  wig = read.table(wig_file_path, skip = 3)
  
  # Combine negative and positive strand counts
  wig_counts = wig[,1:6]        #test2=wig_counts
  wig_counts[,2:6] = wig_counts[,2:6] + wig[,9:13]
  colnames(wig_counts) = c( "Pos", "A", "C", "G", "T", "N" )
  
  #Get total read counts for each position
  tot_r_counts = rowSums(wig_counts[,2:6]) 
  #Delineate the number of SNPs at a position
  b = rowSums(wig_counts[,2:5]>0)
  
  # Check and keep positions with only two alleles
  f =  b == 2
  
  # Calculate ratios and merge data
  d = wig_counts[f,2:5]/tot_r_counts[f]
  d[d==0]=NA 
  d = as.matrix(d)
  
  test3.mask = wig_counts[ f ,1:5]  #Positions and base counts
  test3.mask[test3.mask < 10] = 0       #Set base counts lees than 10 to 0, excluding these.  READ COUNT FILTER##########################################
  b3.mask = rowSums(test3.mask[,2:5]>0) #Filtered, delineate number of SNPS at a position
  
  #If there's no SNPs present, return NA
  if(sum(b3.mask == 2) <= 1){return(NA)}
  
  # Skip polymorphic sites (ie more than SNPs)
  m = match(test3.mask[b3.mask==2,1], vcf[,2])  #Match the positions that have only 2 alleles to their chromosome positions
  f.t = !is.na(m)                               #The ones that matched
  f.v = m[f.t] 
  
  #If there are no heterozygous SNPs with mapped reads on two alleles, skip
  if(length(m) <= 1){return(NA)}
  
  # Tidy up, find ref and alt
  if(dim(vcf[f.v,4:5])[1] <= 1){return(NA)}  #only one heterozygous SNP
  temp.vcf = cbind(d[b3.mask==2,][f.t,], vcf[f.v,4:5]) #get the ratios of positions that only have 2 alleles, and add on the Ref and Alt allele from the VCF file
  
  ref1 = as.character(temp.vcf[,5]) == "A"             #Grab all the reference A's, and so on
  ref2 = as.character(temp.vcf[,5]) == "C"
  ref3 = as.character(temp.vcf[,5]) == "G"
  ref4 = as.character(temp.vcf[,5]) == "T"
  
  alt1 = as.character(temp.vcf[,6]) == "A"             #Grab all the alternative A's, and so on
  alt2 = as.character(temp.vcf[,6]) == "C"
  alt3 = as.character(temp.vcf[,6]) == "G"
  alt4 = as.character(temp.vcf[,6]) == "T"
  
  
  ref.only = c(temp.vcf[ref1,1], temp.vcf[ref2,2], temp.vcf[ref3,3], temp.vcf[ref4,4])  #Combine all the reference SNPs ratios
  alt.only = c(temp.vcf[alt1,1], temp.vcf[alt2,2], temp.vcf[alt3,3], temp.vcf[alt4,4])  #Combine all the alternative SNPs ratios
  
  temp2.vcf = cbind(test3.mask[b3.mask==2,2:5][f.t,], vcf[f.v,4:5])                     #Grab the filtered SNPs and their alleles, for only those with two alleles
  temp3.vcf = temp2.vcf[,1:2]                                                           #Set up structure to hold the combined Ref / Alt counts per position of the filtered SNPS
  temp3.vcf = temp3.vcf*0
  colnames(temp3.vcf) = c('ref_counts','alt_counts')
  
  #Starting filling in the ref and alt counts
  #Produces the ref and alt read counts for the heterozygous SNPs
  temp3.vcf[ref1,1] = temp2.vcf[ref1,1]
  temp3.vcf[ref2,1] = temp2.vcf[ref2,2]
  temp3.vcf[ref3,1] = temp2.vcf[ref3,3]
  temp3.vcf[ref4,1] = temp2.vcf[ref4,4]
  
  temp3.vcf[alt1,2] = temp2.vcf[alt1,1]
  temp3.vcf[alt2,2] = temp2.vcf[alt2,2]
  temp3.vcf[alt3,2] = temp2.vcf[alt3,3]
  temp3.vcf[alt4,2] = temp2.vcf[alt4,4]
  
  #filtering out Variants with , in them, they are multiple alleles missed previously
  allele.f1 = sapply(1:length(temp2.vcf$V4),function(i) str_count(as.character(temp2.vcf$V4)[i], ',') )    #Ref alleles with commas
  allele.f2 = sapply(1:length(temp2.vcf$V5),function(i) str_count(as.character(temp2.vcf$V5)[i], ',') )    #Alt alleles with commas
  allele.f = !(allele.f1 == 1 | allele.f2 == 1)                                                            #Get rid of them
  
  temp2.vcf = temp2.vcf[allele.f, ]
  temp3.vcf = temp3.vcf[allele.f, ]
  
  vcf_2 = cbind(test3.mask[b3.mask==2,1][f.t][allele.f], temp3.vcf, temp2.vcf)
  colnames(vcf_2) = c('Pos','ref_counts','alt_counts','A','C','G','T', 'ref','alt')
  
  #Any remaining SNP with 0 reads for the reference or alternate allele indicates a disagreement between the VCF identified alt and ref allele and the bases with read pileups
  #Just filter them out
  vcf_2 = vcf_2[!(vcf_2$ref_counts == 0 | vcf_2$alt_counts == 0), ]
  
  return(vcf_2)
  
}