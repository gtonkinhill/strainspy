testthat::test_that("strainspy wrapper fails with clear message for missing design variable", {
  d <- ss_test_data()

  bad_meta <- d$meta
  # Deliberately remove a design term to force an informative error
  bad_meta$Age_at_collection <- NULL

  meta_file <- tempfile(fileext = ".csv")
  sylph_file <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")

  utils::write.csv(
    bad_meta,
    meta_file,
    row.names = FALSE
  )

  testthat::expect_error(
    strainspy::strainspy(
      meta_path = meta_file,
      sylph_path = sylph_file,
      design_formula = "~ Case_status + Age_at_collection",
      coef = 2,
      nthreads = 1,
      alpha = 1
    ),
    "Metadata is missing predictor column\\(s\\)"
  )
})

testthat::test_that("strainspy wrapper requires design_formula and coef", {
  d <- ss_test_data()
  meta_file <- tempfile(fileext = ".csv")

  utils::write.csv(d$meta, meta_file, row.names = FALSE)

  testthat::expect_error(
    strainspy::strainspy(
      meta_path = meta_file,
      sylph_path = d$sylph_path,
      coef = 2
    ),
    "design_formula"
  )

  testthat::expect_error(
    strainspy::strainspy(
      meta_path = meta_file,
      sylph_path = d$sylph_path,
      design_formula = "~ Case_status + Age_at_collection"
    ),
    "`coef` is required"
  )
})

testthat::test_that("strainspy wrapper runs with generic configurable workflow", {
  d <- ss_test_data()

  meta <- d$meta
  meta_file <- tempfile(fileext = ".csv")
  out_file <- tempfile(fileext = ".csv")

  utils::write.csv(
    meta,
    meta_file,
    row.names = FALSE
  )

  out <- suppressWarnings(
    strainspy::strainspy(
      meta_path = meta_file,
      sylph_path = system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy"),
      taxonomy_path = system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy"),
      output_path = out_file,
      design_formula = "~ Case_status + Age_at_collection",
      coef = 2,
      nthreads = 1,
      alpha = 1,
      return_vcov = FALSE
    )
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true(file.exists(out_file))
  testthat::expect_true(all(c("coefficient", "p_adjust") %in% colnames(out)))
})
