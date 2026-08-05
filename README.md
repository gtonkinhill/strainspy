
<!-- README.md is generated from README.Rmd. Please only edit the Rmd file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/gtonkinhill/strainspy/workflows/R-CMD-check-hard/badge.svg)](https://github.com/gtonkinhill/strainspy/actions)
[![Databases](https://zenodo.org/badge/DOI/10.5281/zenodo.21796878.svg)](https://zenodo.org/records/21796878)
<!-- badges: end -->

# Strainspy

## Installation

`strainspy` is currently available on github. It can be installed with
`remotes`

``` r
install.packages("remotes")
remotes::install_github("gtonkinhill/strainspy")
```

If you would like to also build the vignette with your installation run:

``` r
remotes::install_github("gtonkinhill/strainspy", build_vignettes = TRUE)
```

## Example

This walkthrough demonstrates a typical `strainspy` analysis and
showcases some of the models and outputs available. Here, we analyse a
200 sample subset of the data described in [Wallen *et al.*
2022](https://doi.org/10.1038/s41467-022-34667-x).

**NOTE:** Be sure to replace the example paths with valid file paths on
your system. The `system.file()` function is used here only for
demonstration purposes in this vignette and will not work outside the
package environment.

``` r
library(strainspy)
library(ggplot2)

example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
example_sylph_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
example_taxonomy_path <- system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy")
```

### Load and filter data

``` r
# Read in metadata
meta_data <- readr::read_csv(example_meta_path)

meta_data$Case_status = factor(meta_data$Case_status)  # Required for visualising
# Read in sylph profile
se <- read_sylph(example_sylph_path, meta_data = meta_data)

# Filter by presence. This will remove any strains that are not present in at
# least 30 samples
se <- filter_by_presence(se, min_nonzero = 30)
#> Retained 472 rows after filtering

# Read in taxonomy
taxonomy <- read_taxonomy(example_taxonomy_path)
```

### Fit the model

``` r
# Create design matrix
design <- as.formula(" ~ Case_status + Age_at_collection")

# Fit a Zero-inflated beta model using the default preset_weak MAP prior
fit <- glmZiBFit(se, design, nthreads = parallel::detectCores())
#> Fitting model... 
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
```

### Summarise and plot the results

``` r
# Get top hits
th = top_hits(fit, coef = 2)
#> Found 2 tophits for Case_statusPD at alpha = 0.05 using holm

# Create Volcano plot
plot_volcano(fit, label = T)
#> Found 472 tophits for Case_statusPD at alpha = 1 using holm
```

<img src="inst/vignette-supp/unnamed-chunk-7-1.png" alt="" width="100%" />

``` r

knitr::kable(th)
```

| Contig_name | Genome_file | coefficient | std_error | p_value | p_adjust | zi_coefficient | zi_std_error | zi_p_value | zi_p_adjust |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| NZ_WSNW01000023.1 Coprobacillus cateniformis strain RCA1-24 10234_10235, whole genome shotgun sequence | GCF_009767585.1 | 0.1667314 | 0.1243133 | 0.1798493 | 1.0000000 | -1.8172219 | 0.4454086 | 0.0000451 | 0.0212661 |
| UREB01000001.1 uncultured Clostridiales bacterium isolate UMGS913 genome assembly, contig: NODE_380_length_41473_cov_4.232942, whole genome shotgun sequence | GCA_900546685.1 | 0.2481419 | 0.0638069 | 0.0001007 | 0.0475203 | -0.3495862 | 0.2874739 | 0.2239611 | 1.0000000 |

*Coprobacillus cateniformis* shows difference in presence (adjusted p =
0.0213). *Clostridiales bacterium* shows difference in identity
(adjusted p \< 0.01).

### Perform post-hoc testing to validate the beta hit

``` r
th = comp_ani_diff_and_posthoc_test(se, fit, th)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

# This NA indicates posthoc testing does not invalidate this beta hit
is.na(th$Comment)
#> [1] TRUE TRUE

th$ANI_Difference
#> [1]        NA 0.2481419
```

In *Clostridiales bacterium*, average ANI (green dashed line) is 0.25%
lower in Controls compared to Parkinson Disease (PD). The distribution
appears different between groups.

## Visualise the distribution of top hits with Case_status

``` r
plot_ani_dist(se, "Case_status", top_hits(fit)$Contig_name)
#> Found 2 tophits for Case_statusPD at alpha = 0.05 using holm
```

<img src="inst/vignette-supp/unnamed-chunk-9-1.png" alt="" width="100%" />

``` r
plot_histogram(se, "Case_status", top_hits(fit)$Contig_name[2], drop_zeros = TRUE,
    fit_spline = T, bins = 54, separate_facets = T)
#> Found 2 tophits for Case_statusPD at alpha = 0.05 using holm
```

<img src="inst/vignette-supp/unnamed-chunk-9-2.png" alt="" width="100%" />

## Incorporate taxonomy

### Generate taxonomy informed Manhattan plot with adjusted p-values

``` r
plot_manhattan(fit, taxonomy = taxonomy)
```

<img src="inst/vignette-supp/unnamed-chunk-10-1.png" alt="" width="100%" />

### Create a traditional Manhattan plot coloured by taxonomy with unadjusted p-values and Bonferroni significance thresholds

``` r
plot_manhattan(fit, taxonomy = taxonomy, aggregate_by_taxa = F)
```

![](inst/vignette-supp/unnamed-chunk-11-1.png)<!-- -->

## Estimate sample richness

This can be used as a proxy for alpha diversity.

``` r
rchns = estimate_sample_richness(se, taxonomy)$species_richness
df_rch = merge(data.frame(run_accession = names(rchns), value = rchns), meta_data,
    by = "run_accession")

ggplot(df_rch, aes(Case_status, value, fill = Case_status)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4) + theme_classic() + ylab("Number of Species")
```

<img src="inst/vignette-supp/unnamed-chunk-12-1.png" alt="" width="100%" />

## Estimate between sample diversity

``` r
divs = strainspy::estimate_beta_diversity(se, phenotype = "Case_status", return_plots = T)

divs$plots$NMDS
```

<img src="inst/vignette-supp/unnamed-chunk-13-1.png" alt="" width="100%" />

## Working with other inputs

### Example using Sourmash output

``` r
example_sourmash_path <- system.file("extdata", "example_sourmash.csv.gz", package = "strainspy")
sm <- read_sourmash(example_sourmash_path, meta_data)
```

All remaining functions are compatible with this output.

**Note:** `strainspy` provides a function to merge `sourmash gather` and
`sourmash search` outputs. See help for details.

### Example using MetaPhlAn output

``` r
example_metaphlan_path <- system.file("extdata", "metaphlan_merged.tsv.gz", package = "strainspy")
example_taxonomy_path <- system.file("extdata", "metaphlan_taxonomy.tsv.gz", package = "strainspy")
mp <- read_metaphlan(example_metaphlan_path, meta_data)
```

All remaining functions are compatible with this output.

**Note:** `strainspy` provides a function to merge `MetaPhlAn` profiles
and generate the taxonomy file. See help for details.

## Citation

To cite strainspy please use
