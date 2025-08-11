
<!-- README.md is generated from README.Rmd. Please only edit the Rmd file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/gtonkinhill/strainspy/workflows/R-CMD-check-hard/badge.svg)](https://github.com/gtonkinhill/strainspy/actions)
<!-- [![DOI](https://zenodo.org/badge/XXXX.svg)](https://zenodo.org/badge/latestdoi/XXXX) -->
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
```

### Fit the model

``` r
# Create design matrix
design <- as.formula(" ~ Case_status + Age_at_collection")

# Fit a Zero-inflated beta model using the default preset_weak MAP prior
fit <- glmZiBFit(se, design, nthreads = parallel::detectCores())
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
```

### Summarise and plot the results

``` r
# Get top hits
top_hits(fit, coef = 2)
#> # A tibble: 3 × 10
#>   Contig_name  Genome_file coefficient std_error p_value p_adjust zi_coefficient
#>   <chr>        <chr>             <dbl>     <dbl>   <dbl>    <dbl>          <dbl>
#> 1 NZ_WSNW0100… GCF_009767…       0.168    0.124  1.75e-1   1              -1.80 
#> 2 UREB0100000… GCA_900546…       0.248    0.0638 9.86e-5   0.0466         -0.349
#> 3 DVOB0100001… GCA_018713…       0.232    0.0598 1.05e-4   0.0497         -0.663
#> # ℹ 3 more variables: zi_std_error <dbl>, zi_p_value <dbl>, zi_p_adjust <dbl>

# Create Volcano plot
plot_volcano(fit, label = T)
```

<img src="inst/vignette-supp/unnamed-chunk-7-1.png" width="100%" />

## Visualise the distribution of top hits with Case_status

### *Coprobacillus cateniformis* shows difference in presence

``` r
plot_ani_dist(se, "Case_status", top_hits(fit)$Contig_name[1], show_points = T)
```

<img src="inst/vignette-supp/unnamed-chunk-8-1.png" width="100%" />

### *Clostridiales bacterium* and *Candidatus Copromorpha excrementipullorum* shows difference in identity

``` r
plot_ani_dist(se, "Case_status", top_hits(fit)$Contig_name[2:3], show_points = T,
    drop_zeros = T, plot_type = "box")
```

<img src="inst/vignette-supp/unnamed-chunk-9-1.png" width="100%" />

## Incorporate taxonomy

### Generate taxonomy informed Manhattan plot with adjusted p-values

``` r
# Read in taxonomy
taxonomy <- read_taxonomy(example_taxonomy_path)

plot_manhattan(fit, taxonomy = taxonomy)
```

<img src="inst/vignette-supp/unnamed-chunk-10-1.png" width="100%" />

### Create a traditional Manhattan plot coloured by taxonomy with unadjusted p-values and Bonferroni significance thresholds

``` r
plot_manhattan(fit, taxonomy = taxonomy, aggregate_by_taxa = F)
```

![](inst/vignette-supp/unnamed-chunk-11-1.png)<!-- -->

## Example using Sourmash output

``` r
example_sourmash_path <- system.file("extdata", "example_sourmash.csv.gz", package = "strainspy")
sm <- read_sourmash(example_sourmash_path, meta_data)
```

All remaining functions are compatible with this output.

**Note:** `strainspy` provides a function to merge `sourmash gather` and
`sourmash search` outputs. See help for details.

## Example using MetaPhlAn output

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
