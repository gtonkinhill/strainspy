
<!-- README.md is generated from README.Rmd. Please only edit the Rmd file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/gtonkinhill/strainspy/workflows/R-CMD-check-hard/badge.svg)](https://github.com/gtonkinhill/strainspy/actions)
<!-- [![DOI](https://zenodo.org/badge/XXXX.svg)](https://zenodo.org/badge/latestdoi/XXXX) -->
<!-- badges: end -->

# Strainspy

## Project Status: Active Development

StrainSpy is currently under active development. APIs and functionality
may change without notice, including backward-incompatible updates.

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

<img src="inst/vignette-supp/unnamed-chunk-7-1.png" width="100%" />

``` r

th
#> # A tibble: 2 × 10
#>   Contig_name  Genome_file coefficient std_error p_value p_adjust zi_coefficient
#>   <chr>        <chr>             <dbl>     <dbl>   <dbl>    <dbl>          <dbl>
#> 1 NZ_WSNW0100… GCF_009767…       0.167    0.124  1.80e-1   1              -1.82 
#> 2 UREB0100000… GCA_900546…       0.248    0.0638 1.01e-4   0.0475         -0.350
#> # ℹ 3 more variables: zi_std_error <dbl>, zi_p_value <dbl>, zi_p_adjust <dbl>
```

*Coprobacillus cateniformis* shows difference in presence (adjusted p =
0.0213). *Clostridiales bacterium* shows difference in identity
(adjusted p \< 0.01).

### Perform post-hoc testing to validate the beta hit

``` r
th = comp_ani_diff_and_posthoc_test(se, fit, th)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

# This NA indicates posthoc testing does not invalidate this beta hit
is.na(th$Comment)
#> [1] TRUE

th$ANI_Difference
#> [1] -0.2481419
```

In *Clostridiales bacterium*, average ANI is 0.25% lower in Controls
compared to Parkinson Disease (PD).

## Visualise the distribution of top hits with Case_status

``` r
plot_ani_dist(se, "Case_status", top_hits(fit)$Contig_name)
#> Found 2 tophits for Case_statusPD at alpha = 0.05 using holm
```

<img src="inst/vignette-supp/unnamed-chunk-9-1.png" width="100%" />

## Incorporate taxonomy

### Generate taxonomy informed Manhattan plot with adjusted p-values

``` r
plot_manhattan(fit, taxonomy = taxonomy)
```

<img src="inst/vignette-supp/unnamed-chunk-10-1.png" width="100%" />

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

<img src="inst/vignette-supp/unnamed-chunk-12-1.png" width="100%" />

## Estimate between sample diversity

``` r
divs = strainspy::estimate_beta_diversity(se, phenotype = "Case_status", return_plots = T)
#> Run 0 stress 0.2209113 
#> Run 1 stress 0.2208914 
#> ... New best solution
#> ... Procrustes: rmse 0.001283966  max resid 0.01049977 
#> Run 2 stress 0.2221602 
#> Run 3 stress 0.2215367 
#> Run 4 stress 0.2224751 
#> Run 5 stress 0.2208967 
#> ... Procrustes: rmse 0.001241622  max resid 0.01047751 
#> Run 6 stress 0.2235561 
#> Run 7 stress 0.2303237 
#> Run 8 stress 0.2209092 
#> ... Procrustes: rmse 0.001346247  max resid 0.01265921 
#> Run 9 stress 0.2213041 
#> ... Procrustes: rmse 0.006201783  max resid 0.08068013 
#> Run 10 stress 0.2241569 
#> Run 11 stress 0.2356318 
#> Run 12 stress 0.220944 
#> ... Procrustes: rmse 0.001679674  max resid 0.01045541 
#> Run 13 stress 0.2209528 
#> ... Procrustes: rmse 0.002181134  max resid 0.0138078 
#> Run 14 stress 0.2215357 
#> Run 15 stress 0.2355633 
#> Run 16 stress 0.2305136 
#> Run 17 stress 0.2232026 
#> Run 18 stress 0.221555 
#> Run 19 stress 0.2231904 
#> Run 20 stress 0.2652676 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>      1: no. of iterations >= maxit
#>     18: stress ratio > sratmax
#>      1: scale factor of the gradient < sfgrmin

divs$plots$NMDS
```

<img src="inst/vignette-supp/unnamed-chunk-13-1.png" width="100%" />

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
