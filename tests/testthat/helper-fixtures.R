ss_test_data <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
    sylph_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
    taxonomy_path <- system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy")
    metaphlan_path <- system.file("extdata", "metaphlan_merged.tsv.gz", package = "strainspy")
    sourmash_path <- system.file("extdata", "example_sourmash.csv.gz", package = "strainspy")

    stopifnot(
      nzchar(meta_path),
      nzchar(sylph_path),
      nzchar(taxonomy_path),
      nzchar(metaphlan_path),
      nzchar(sourmash_path)
    )

    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    meta$Case_status <- as.factor(meta$Case_status)

    se <- strainspy::read_sylph(sylph_path, meta_data = meta)
    se_filtered <- strainspy::filter_by_presence(se, min_nonzero = 30)
    taxonomy <- strainspy::read_taxonomy(taxonomy_path)

    cache <<- list(
      meta = meta,
      meta_path = meta_path,
      sylph_path = sylph_path,
      taxonomy_path = taxonomy_path,
      metaphlan_path = metaphlan_path,
      sourmash_path = sourmash_path,
      se = se,
      se_filtered = se_filtered,
      taxonomy = taxonomy
    )

    cache
  }
})

ss_fit_zib <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- suppressWarnings(
      strainspy::glmZiBFit(
        d$se_filtered[1:20, ],
        design,
        nthreads = 1,
        return_vcov = FALSE
      )
    )

    cache
  }
})

ss_fit_zib_readme <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- suppressWarnings(
      strainspy::glmZiBFit(
        d$se_filtered,
        design,
        nthreads = 1,
        return_vcov = FALSE
      )
    )

    cache
  }
})

ss_fit_qb <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- strainspy::glmQBFit(
      d$se_filtered[1:30, ],
      design,
      nthreads = 1
    )

    cache
  }
})
