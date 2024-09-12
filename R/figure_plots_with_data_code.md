Werner, Hover, and Gillis. Population variability in X-chromosome
inactivation across 10 mammalian species. 2024. Code and data for
figures.
================
Jonathan Werner
2024-09-12

This markdown file contains the code for generating all plots for the
publication: Werner, Hover, and Gillis. Population variability in
X-chromosome inactivation across 10 mammalian species. 2024.

All data for these plots is provided at:
<https://github.com/gillislab/cross_mammal_xci/tree/main/R/data_for_plots>

``` r
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(MetBrewer)
species_colors = met.brewer("Tiepolo", n=10,type="continuous")

#Going by the ranking of estimated cell counts
species_palette = c('Macaca' = species_colors[1], 'Rat' = species_colors[2], 'Pig' = species_colors[3],
                    'Goat' = species_colors[4], 'Horse' = species_colors[5], 'Sheep' = species_colors[6],
                    'Human' = species_colors[8], 'Mouse' = species_colors[7], 'Cow' = species_colors[9], 'Dog' = species_colors[10])

species_palette
```

    ##    Macaca       Rat       Pig      Goat     Horse     Sheep     Human     Mouse 
    ## "#802417" "#B1572F" "#C77F3D" "#D69F4D" "#D9B05B" "#72763F" "#3B7074" "#3D5F49" 
    ##       Cow       Dog 
    ## "#437E96" "#17486F"

## Figure 1C

``` r
#Contains the species_sample_df dataframe
load('/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/species_sample_df.Rdata')


filter_label = c(rep('Annotated female', length = 10), rep('>= 10 well-powered SNPs', length = 10))
ggplot(species_sample_df, aes(x = species, y = sample_number, alpha = filter_label, fill = species)) + geom_bar(position="dodge", stat="identity") + coord_flip() +
  ylab('# of bulk RNA-seq samples') +
  scale_fill_manual(values = species_palette, guide = guide_legend(reverse = TRUE)) + scale_alpha_manual(values = c('Annotated female' = 1, '>= 10 well-powered SNPs' = .5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = 'white'), panel.border = element_rect(color = "black", fill=NA),  
      axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
      axis.title.y = element_text(size=12), axis.title.x = element_text(size=12)) 
```

![](figure_plots_with_data_code_files/figure-gfm/species_sample-1.png)<!-- -->
