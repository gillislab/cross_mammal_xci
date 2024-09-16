Description of R notebooks

get_exp_vcf.R - R script for combining the .vcf and .wig output from the snakemake pipeline to produce tables of variants with their allele-specific read counts

cross_mammal_data.Rmd - R notebook compiling sample metadata from SRA

exploring_SPECIES_samples.Rmd - R notebook per species for compiling the X-chromosome .vcfs, filtering and annotating SNPs, filtering X-escape, and sample-level XCI ratio modeling. Plots for Supplemental Figures 1, 2, 4

exploring_SPECIES_samples_chrA_and_B.Rmd - R noteboot per species for compiling the autosome .vcfs, filtering and annotating SNPs, and sample-level autosome ratio modeling. Analysis for Supplemental Figure 5.

aggregate_gtex_data.Rmd - R notebook that aggregates the tissue vcfs per donor for the human GTEx data

autosome_analysis.Rmd - R notebook comparing autosomal and X-chromosome imbalances and adding a filter for samples with global allelic imbalances. Plots for Supplemental Figure 5.

crossSpecies_cell_num_modeling.Rmd - R notebook that contains analysis for the X-linked heterozygosity and population XCI ratio modeling results. Plots for Figures 2-3 and Supplemental Figures 6-7.

xci_predict_variant_analysis.Rmd - R notebook that contains analysis for quantifying the association between individual X-linked variants and extreme XCI ratios across species. Plots for Figure 4 and Supplemental Figures 8, 9, and 10.

human_tissue_specific_XCI_variant.Rmd - R notebook containing analysis for the human XCI ratio-variant associations, Supp. Fig. 9 A

fig_1_schematic_plots.Rmd - R notebook for miscellaneous plots used in the Figure 1 schematic

mouse_sex_classification.Rmd - R notebook containing analysis for annotating the sex of the mouse samples, comparing aggregate chrY and chrX counts per sample, Supp. Figure 3.
