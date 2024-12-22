#' Read Sylph query or profile output file and create a SummarizedExperiment object
#'
#' This function reads a Sylph output file (either query or profile), parses the data, and returns a SummarizedExperiment object.
#' Optionally, metadata can be provided, which will be loaded into the `colData` of the SummarizedExperiment object.
#'
#' @param file_path Character. Path to the Sylph query or profile output file (tab-separated format).
#' @param meta_data data.frame. A tibble or data frame containing sample metadata.
#' The **first column must contain sample names** that match exactly with the sample names in the Sylph output.
#' @param clean_names Logical. If `TRUE`, file paths will be stripped of their directory path and file extension,
#' leaving only the base file name. Defaults to `TRUE`.
#'
#' @return A SummarizedExperiment object with the following components:
#'   \item{assays}{A matrix containing numeric features such as `Adjusted_ANI`, `Taxonomic_abundance`, `Median_cov`, etc.}
#'   \item{rowData}{DataFrame containing metadata for each sequence, such as `Contig_name` and `Genome_file`.}
#'   \item{colData}{DataFrame containing information for each sample (derived from the `Sample_file` field).}
#'
#' @importFrom readr read_tsv
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Read a Sylph file (query or profile) into a SummarizedExperiment object
#'   example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainseekr")
#'   se <- read_sylph(example_path)
#'   # View the SummarizedExperiment
#'   se
#'   # View the assays (numerical matrix)
#'   assay(se)
#'   # View the rowData (metadata for contigs)
#'   rowData(se)
#'   # View the colData (metadata for samples)
#'   colData(se)
#' }
#' \dontrun{
#'   # Read a Sylph file (query or profile) with associated metadata into a SummarizedExperiment object
#'   example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainseekr")
#'   example_meta <- readr::read_csv(example_meta_path)
#'   example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainseekr")
#'   se <- read_sylph(example_path, example_meta)
#' }
read_sylph <- function(file_path, meta_data=NULL, clean_names = TRUE) {

  # Check input argument
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("`file_path` must be a character string specifying the path to the Sylph output file.")
  }

  if (!file.exists(file_path)) {
    stop(paste0("The file '", file_path, "' does not exist. Please provide a valid file path."))
  }

  # Read the Sylph output file using readr::read_tsv
  sylph_data <- readr::read_tsv(
    file_path,
    na = c("", "NA"),
    col_types = readr::cols(
      Sample_file = readr::col_character(),
      Genome_file = readr::col_character(),
      Taxonomic_abundance = readr::col_double(),
      Sequence_abundance = readr::col_double(),
      Adjusted_ANI = readr::col_double(),
      True_cov = readr::col_double(),
      `ANI_5-95_percentile` = readr::col_character(),
      Eff_lambda = readr::col_character(),
      `Lambda_5-95_percentile` = readr::col_character(),
      Median_cov = readr::col_double(),
      Mean_cov_geq1 = readr::col_double(),
      Containment_ind = readr::col_character(),
      Naive_ANI = readr::col_double(),
      kmers_reassigned = readr::col_double(),
      Contig_name = readr::col_character()
    )
  )

  # Identify if the file is a query or profile output based on the presence of specific columns
  is_profile_output <- all(c(
    "Taxonomic_abundance",
    "Sequence_abundance",
    "True_cov",
    "kmers_reassigned"
  ) %in% colnames(sylph_data))

  if (is_profile_output) {
    message("Detected Sylph profile output file.")
    required_columns <- c(
      "Sample_file", "Genome_file", "Taxonomic_abundance", "Sequence_abundance",
      "Adjusted_ANI", "True_cov", "ANI_5-95_percentile",
      "Eff_lambda", "Lambda_5-95_percentile", "Median_cov",
      "Mean_cov_geq1", "Containment_ind", "Naive_ANI",
      "kmers_reassigned", "Contig_name"
    )
  } else {
    message("Detected Sylph query output file.")
    required_columns <- c(
      "Sample_file", "Genome_file", "Adjusted_ANI", "Eff_lambda",
      "ANI_5-95_percentile", "Lambda_5-95_percentile", "Median_cov",
      "Mean_cov_geq1", "Containment_ind", "Naive_ANI", "Contig_name"
    )
  }

  # Validate that the data contains the expected columns
  missing_columns <- setdiff(required_columns, colnames(sylph_data))
  if (length(missing_columns) > 0) {
    stop(paste0("The following required columns are missing from the Sylph file: ",
                paste(missing_columns, collapse = ", ")))
  }

  # Deal with file names and extensions
  if (clean_names){
    sylph_data$Sample_file <- tools::file_path_sans_ext(basename(sylph_data$Sample_file),
                                                        compression = TRUE)
    sylph_data$Genome_file <- tools::file_path_sans_ext(basename(sylph_data$Genome_file),
                                                        compression = TRUE)
  }

  # Merge metadata if provided
  if (is.null(meta_data)){
    # Generate colData
    col_data <- as.matrix(unique(sylph_data$Sample_file))
    rownames(col_data) <- col_data[,1]

  } else {
    col_data <- data.frame(Sample_file = unique(sylph_data$Sample_file))

    # Extract unique sample names
    sylph_samples <- unique(sylph_data$Sample_file)
    meta_samples <- unique(meta_data[[1]])

    # Warn about mismatched samples
    if (length(missing_from_meta <- setdiff(sylph_samples, meta_samples)) > 0) {
      stop("The following samples from 'sylph_data' are not in 'meta_data': ", paste(missing_from_meta, collapse = ", "))
    }

    if (length(missing_from_sylph <- setdiff(meta_samples, sylph_samples)) > 0) {
      stop("The following samples from 'meta_data' are not in 'sylph_data': ", paste(missing_from_sylph, collapse = ", "))
    }

    col_data <- base::merge(col_data, meta_data,
                               by.x = "Sample_file",
                               by.y = names(meta_data)[1],
                               all.x = TRUE) # Keeps all Sample_file entries from sylph_data

    # Generate colData
    rnames <- col_data[,1]
    rownames(col_data) <- rnames

  }

  # Extract row metadata (rowData)
  row_data <- as.matrix(unique(sylph_data[, c(
    "Contig_name", "Genome_file"
  )]))
  rownames(row_data) <- row_data[,1]

  # Extract numeric data to create the "assay" matrix
  if (is_profile_output) {
    numeric_columns <- c(
      "Taxonomic_abundance", "Sequence_abundance",
      "Adjusted_ANI", "Eff_lambda", "Median_cov"
    )
  } else {
    numeric_columns <- c("Adjusted_ANI", "Eff_lambda", "Median_cov")
  }

  assay_data <- lapply(numeric_columns, function(var){
    m <- matrix(0, nrow = nrow(row_data), ncol = nrow(col_data),
           dimnames = list(rownames(row_data), rownames(col_data)))
    m[cbind(sylph_data$Contig_name, sylph_data$Sample_file)] <- sylph_data[[var]]
    return(m)
  })

  # Create the SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = setNames(assay_data, numeric_columns),
    rowData = row_data,
    colData = col_data
  )

  # Return the SummarizedExperiment object
  return(se)
}
