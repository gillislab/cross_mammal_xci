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
library(ggrepel)
```

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

## Supp. Figure 1B

``` r
load('/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/num_snps_df.Rdata')

ggplot(num_snps_df, aes(x = species, y = snp_counts, fill = species)) + geom_violin(scale = 'width') + scale_fill_manual(values = species_palette) +
  ylab('Number of well-powered SNPs per sample')
```

![](figure_plots_with_data_code_files/figure-gfm/num_snps-1.png)<!-- -->

## Supp. Figure 6A

``` r
#contains the species_squared_error_df dataframe
load(file = '/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/species_squared_error_df')


plot(species_squared_error_df$n_vec, species_squared_error_df$horse_error, main = 'Horse squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$horse_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-1.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$dog_error, main = 'dog squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$dog_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-2.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$cow_error, main = 'cow squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$cow_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-3.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$sheep_error, main = 'sheep squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$sheep_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-4.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$goat_error, main = 'goat squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$goat_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-5.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$rat_error, main = 'rat squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$rat_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-6.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$mouse_error, main = 'mouse squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$mouse_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-7.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$macaca_error, main = 'macaca squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$macaca_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-8.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$pig_error, main = 'pig squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$pig_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-9.png)<!-- -->

``` r
plot(species_squared_error_df$n_vec, species_squared_error_df$human_error, main = 'human squared error', xlab = 'Cell number estimate', ylab = 'Sum squared error')
abline(v = species_squared_error_df$n_vec[which.min(species_squared_error_df$human_error)], col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_error_dot_plots-10.png)<!-- -->

## Figure 2A

``` r
#Contains the horse_skew_and_stats_df, dog_skew_and_stats_df, cow_skew_and_stats_df, sheep_skew_and_stats_df, goat_skew_and_stats_df, rat_skew_and_stats_df, 
#mouse_skew_and_stats_df, macaca_skew_and_stats_df, pig_skew_and_stats_df, agg_gtex_skew_and_stats_df dataframes
load(file = '/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/all_species_skew_and_stats_df.Rdata')

all_skew_and_stats = list(horse_skew_and_stats_df, dog_skew_and_stats_df, cow_skew_and_stats_df, sheep_skew_and_stats_df, 
                          goat_skew_and_stats_df, rat_skew_and_stats_df, mouse_skew_and_stats_df, macaca_skew_and_stats_df, pig_skew_and_stats_df, agg_gtex_skew_and_stats_df)
 
names(all_skew_and_stats) = c('Horse','Dog','Cow','Sheep','Goat','Rat','Mouse','Macaca','Pig', 'Human')


# Functions to fold and unfold ratios 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 


p = .5
num_bins = 100
binned_x_axis = seq(0,1,1/num_bins)


for(j in 1:length(all_skew_and_stats)){

  fitted_normal_n = species_squared_error_df$n_vec[which.min(species_squared_error_df[ ,j+1])] #first column is the n size vector, skip that
  sigma = sqrt(p*(1-p)/fitted_normal_n)
  dnorm_x = 0:1000/1000
  
  fitted_normal = dnorm(dnorm_x, mean=p, sd = sigma)
  
  theoretical_norms = list()
  n = c(8, 32)
  for(i in 1:length(n)){ 
    sigma_t = sqrt(p*(1-p)/n[i])
    theoretical_norms[[i]] = sigma_t
  }
  
  species_name = names(all_skew_and_stats)[j]
  species_skew_and_stat_df = all_skew_and_stats[[j]]
  specie_color = species_palette[names(species_palette)== species_name]
  
  #Greying out the .4-.6 skews
  fill_vec = rep(specie_color, length(binned_x_axis))
  fill_vec[binned_x_axis >= 0.4 & binned_x_axis <= 0.6] = 'grey'
  color_vec = rep(specie_color, length(binned_x_axis))
  color_vec[binned_x_axis >= 0.4 & binned_x_axis <= 0.6] = 'grey'
  
  
  if(species_name == 'Human'){
    skews = unfold(agg_gtex_skew_and_stats_df$est_skew[agg_gtex_skew_and_stats_df$num_good_snps >= 10])
    df = data.frame(skews = skews)
  }else{
    skews = unfold(filter(species_skew_and_stat_df, num_good_snps_no_escape >= 10 & autosome_imbalance == F)$est_skew_no_escape)
    df = data.frame(skews = skews)
  }
  
  g = ggplot(df, aes(x = skews)) + geom_histogram(aes(y = ..density..), fill = fill_vec, color = color_vec, binwidth = 1/num_bins) +
    xlim(0,1)+ xlab(sprintf('%s XCI ratios', species_name)) + ggtitle(sprintf('%s population XCI ratio distribution', species_name)) +
    stat_function(fun = dnorm, args = list(mean = .5, sd = theoretical_norms[[1]]), n = num_bins, size = .5, alpha = .5, aes(linetype = 'dashed')) +
    stat_function(fun = dnorm, args = list(mean = .5, sd = theoretical_norms[[2]]), n = num_bins, size = .5, alpha = .5, aes(linetype = 'dotted')) + 
    stat_function(fun = dnorm, args = list(mean = .5, sd = sigma), n = num_bins, size = .75, aes(linetype = 'solid')) +
    scale_linetype_manual(name="Estimated cell number",values=c('dashed','dotted','solid'),
                          labels = c('8 cells','32 cells',sprintf('Estimated %i cells', fitted_normal_n)),
                          breaks=c('dashed','dotted','solid'), 
                          guide = guide_legend(override.aes = list(size = c(.5, .5, 2)))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = 'white'), panel.border = element_rect(color = "black", fill=NA, size=1),  
          axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), 
          legend.position = c(.8,.75))  
  print(g)
  
}
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## ℹ Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The dot-dot notation (`..density..`) was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `after_stat(density)` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-1.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-2.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-3.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-4.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-5.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-6.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-7.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-8.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-9.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](figure_plots_with_data_code_files/figure-gfm/XCI_population_distributions-10.png)<!-- -->

## Figure 2B

``` r
#Contains the all_ci_df dataframe
load('/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/all_ci_df.Rdata')
#Species ordered by evolutionary relationships
all_ci_df$species = factor(all_ci_df$species, levels = c('Dog','Horse','Cow','Goat','Sheep','Pig','Rat','Mouse','Human','Macaca'))
ggplot(all_ci_df, aes(y = species, x = log2(cell_num_est), color = species)) + geom_point(size = 3) + geom_errorbar(aes(xmin = log2(lower_ci), xmax = log2(upper_ci))) +
  scale_color_manual(values = species_palette) + xlim(2, 7) + xlab('Cell divisions (log2 cell # estimate)') + ggtitle('Chromosome X') + ylab('') +
  scale_y_discrete(limits=rev)
```

![](figure_plots_with_data_code_files/figure-gfm/cell_num_est_confidence-1.png)<!-- -->

## Supp. Figure 6 B and C

``` r
#Contains the cell_num_stats_df dataframe
load( file = '/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/cell_num_stats_df.Rdata')

ori_corr = cor(cell_num_stats_df$estimated_cell_num, cell_num_stats_df$sample_number, method = 'spearman')
ggplot(cell_num_stats_df, aes(y = log10(sample_number), x = estimated_cell_num, color = species)) + geom_point(size = 3) +
  scale_color_manual(values = species_palette) + ylab('log10( Total sample number )') + xlab('Cell # estimate') + ggtitle(sprintf('Correlation: %1.3f', ori_corr))
```

![](figure_plots_with_data_code_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ori_corr = cor(cell_num_stats_df$estimated_cell_num, cell_num_stats_df$mean_snp, method = 'spearman')
ggplot(cell_num_stats_df, aes(y =  mean_snp, x = estimated_cell_num, color = species)) + geom_point(size = 3) +
  scale_color_manual(values = species_palette)  + ylab('Mean SNP # per sample') + xlab('Cell # estimate') + ggtitle(sprintf('Correlation: %1.3f', ori_corr))
```

![](figure_plots_with_data_code_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

## Supp. Figure 6 D

``` r
#Contains the compare_ci_df dataframe
load( file = '/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/compare_ci_df.Rdata')

ggplot(compare_ci_df, aes(y = species, x = log2(cell_num_est), color = color_label, alpha = label)) + geom_point(size = 3) + 
  geom_errorbar(aes(xmin = log2(lower_ci), xmax = log2(upper_ci))) +
  scale_alpha_manual(values = c('Downsampled estimate' = .5, 'Original estimate' = 1)) +
  scale_color_manual(values = species_palette) + xlim(2, 7) + xlab('Cell divisions (log2 cell # estimate)') + ggtitle('Chromosome X') + ylab('') +
  scale_y_discrete(limits=rev)
```

![](figure_plots_with_data_code_files/figure-gfm/downsampled_cell_num_ests-1.png)<!-- -->

## Figure 3C and Supp. Figure 7 A and B

``` r
#Contains the cow_het_df, rat_het_df, mouse_het_df, macaca_het_df, pig_het_df, sheep_het_df, goat_het_df, dog_het_df,horse_het_df, human_het_df dataframes
load( '/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/all_species_het_dfs.Rdata') 

all_het_dfs = list(cow_het_df, rat_het_df, mouse_het_df, macaca_het_df, pig_het_df, sheep_het_df, goat_het_df, dog_het_df,horse_het_df, human_het_df)
all_species = c('cow','rat','mouse','macaca','pig','sheep','goat','dog','horse', 'human')

for(i in 1:length(all_het_dfs)){

  het_df = all_het_dfs[[i]]
  species_name = all_species[i]
  corr = cor(het_df$het_fraction, het_df$est_skew_no_escape, method = 'spearman', use = 'complete.obs')
  p1 = ggplot(het_df, aes(x = log10(het_fraction), y = est_skew_no_escape)) + geom_bin2d(bins = 25) + scale_fill_continuous(type = "viridis") +
    ylim(.5, 1) + ylab('Estimated XCI ratio') + xlab('Heterozygousity') + ggtitle(sprintf('%s corr: %0.3f', species_name, corr))
  print(p1)
  
  corr = cor(het_df$het_fraction, het_df$est_var_no_escape, method = 'spearman', use = 'complete.obs')
  p2 = ggplot(het_df, aes(x = log10(het_fraction), y = est_var_no_escape)) + geom_bin2d(bins = 25) + scale_fill_continuous(type = "viridis") +
    ylab('Estimated XCI SD') + xlab('Heterozygousity') + ggtitle(sprintf('%s corr: %0.3f',species_name, corr))
  print(p2)

}
```

    ## Warning: Removed 595 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

    ## Warning: Removed 21 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-1.png)<!-- -->

    ## Warning: Removed 595 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-2.png)<!-- -->

    ## Warning: Removed 303 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).
    ## Removed 21 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-3.png)<!-- -->

    ## Warning: Removed 303 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-4.png)<!-- -->

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-5.png)<!-- -->![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-6.png)<!-- -->

    ## Warning: Removed 232 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

    ## Warning: Removed 20 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-7.png)<!-- -->

    ## Warning: Removed 232 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-8.png)<!-- -->

    ## Warning: Removed 55 rows containing non-finite outside the scale range (`stat_bin2d()`).
    ## Removed 20 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-9.png)<!-- -->

    ## Warning: Removed 55 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-10.png)<!-- -->

    ## Warning: Removed 59 rows containing non-finite outside the scale range (`stat_bin2d()`).
    ## Removed 20 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-11.png)<!-- -->

    ## Warning: Removed 59 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-12.png)<!-- -->

    ## Warning: Removed 39 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-13.png)<!-- -->

    ## Warning: Removed 39 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-14.png)<!-- -->

    ## Warning: Removed 7 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-15.png)<!-- -->

    ## Warning: Removed 7 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-16.png)<!-- -->

    ## Warning: Removed 8 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-17.png)<!-- -->

    ## Warning: Removed 8 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-18.png)<!-- -->

    ## Warning: Removed 342 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

    ## Warning: Removed 23 rows containing missing values or values outside the scale range
    ## (`geom_tile()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-19.png)<!-- -->

    ## Warning: Removed 342 rows containing non-finite outside the scale range
    ## (`stat_bin2d()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_and_XCI-20.png)<!-- -->

## Figure 3A

``` r
species_het_fracs = c(cow_het_df$het_fraction, rat_het_df$het_fraction, mouse_het_df$het_fraction, macaca_het_df$het_fraction, 
                      pig_het_df$het_fraction, sheep_het_df$het_fraction, goat_het_df$het_fraction, dog_het_df$het_fraction,
                      horse_het_df$het_fraction, human_het_df$het_fraction)

species_labels = c(rep('Cow', length = nrow(cow_het_df)), rep('Rat', length = nrow(rat_het_df)), rep('Mouse', length = nrow(mouse_het_df)), 
                   rep('Macaca', length = nrow(macaca_het_df)), 
                   rep('Pig', length = nrow(pig_het_df)), rep('Sheep', length = nrow(sheep_het_df)), rep('Goat', length = nrow(goat_het_df)), 
                   rep('Dog', length = nrow(dog_het_df)), rep('Horse', length = nrow(horse_het_df)),  rep('Human', length = nrow(human_het_df)))

all_species_het_frac_df = data.frame(het_fracs = species_het_fracs, label = species_labels)

median_order_df = all_species_het_frac_df %>% group_by(label) %>% summarize(median = median(het_fracs, na.rm = T)) %>% arrange(desc(median))
all_species_het_frac_df$label = factor( all_species_het_frac_df$label, levels = c(median_order_df$label))

ggplot(all_species_het_frac_df, aes(x = label, y = log10(het_fracs), fill = label)) + geom_violin(scale = 'width') + geom_boxplot(outlier.shape = NA, width = .5, show.legend = F ) +
  ylab('Sample heterozygousity') + xlab('Species') +
  scale_fill_manual(values = species_palette, breaks = levels(all_species_het_frac_df$label))
```

    ## Warning: Removed 1395 rows containing non-finite outside the scale range
    ## (`stat_ydensity()`).

    ## Warning: Removed 1395 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_violins-1.png)<!-- -->

## Figure 3B

``` r
load( file = '/home/werner/projects/cross_species_XCI/final_plots/R/data_for_plots/het_corr_df_v2.Rdata')

ggplot(het_corr_df, aes(x = corrs_label, y = corrs)) + geom_violin(scale = 'width', fill = 'grey25', alpha = .5) + geom_point(aes(color = species_label), size = 5) +
  scale_color_manual(values = species_palette, breaks = rev(levels(all_species_het_frac_df$label))) +
  scale_fill_manual(values = species_palette, breaks = rev(levels(all_species_het_frac_df$label))) +
  geom_label_repel(aes(label = species_label, fill = species_label), color = 'white',box.padding = 0.5) +
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  ylab('Correlation with sample heterozygousity (spearman)') + xlab('')
```

![](figure_plots_with_data_code_files/figure-gfm/heterozygosity_corrs_with_xci-1.png)<!-- -->