testthat::test_that("example readers return valid SummarizedExperiment objects", {
	d <- ss_test_data()

	testthat::expect_s4_class(d$se, "SummarizedExperiment")
	testthat::expect_gt(nrow(d$se), 0)
	testthat::expect_gt(ncol(d$se), 0)
	testthat::expect_true("Case_status" %in% colnames(SummarizedExperiment::colData(d$se)))

	mp <- strainspy::read_metaphlan(d$metaphlan_path, meta_data = d$meta)
	sm <- strainspy::read_sourmash(d$sourmash_path, meta_data = d$meta)

	testthat::expect_s4_class(mp, "SummarizedExperiment")
	testthat::expect_s4_class(sm, "SummarizedExperiment")
	testthat::expect_gt(nrow(mp), 0)
	testthat::expect_gt(nrow(sm), 0)
})

testthat::test_that("reader input validation catches missing files", {
	bad_path <- tempfile(fileext = ".tsv.gz")

	testthat::expect_error(strainspy::read_sylph(bad_path), "does not exist")
	testthat::expect_error(strainspy::read_metaphlan(bad_path), "does not exist")
	testthat::expect_error(strainspy::read_sourmash(bad_path), "does not exist")
})

testthat::test_that("filter_by_presence enforces argument checks and reduces rows", {
	d <- ss_test_data()

	testthat::expect_error(
		strainspy::filter_by_presence(d$se, min_nonzero = 1.5),
		"`min_zero` must be an integer"
	)

	filtered <- strainspy::filter_by_presence(d$se, min_nonzero = 50)
	testthat::expect_s4_class(filtered, "SummarizedExperiment")
	testthat::expect_lte(nrow(filtered), nrow(d$se))
})

testthat::test_that("extract_strainspy returns expected long format columns", {
	d <- ss_test_data()
	contigs <- rownames(d$se_filtered)[1:3]

	out <- strainspy::extract_strainspy(
		d$se_filtered,
		variables = c("Case_status", "Age_at_collection"),
		contigs = contigs,
		taxonomy = d$taxonomy
	)

	testthat::expect_s3_class(out, "data.frame")
	testthat::expect_true(all(c("sample", "Case_status", "Age_at_collection", "Contig_name", "ANI") %in% colnames(out)))
	testthat::expect_true(all(unique(out$Contig_name) %in% contigs))
})