#' Read Sylph query or profile output file and create a SummarizedExperiment object
#'
#' This function reads a Sylph output file (either query or profile), parses the data, and returns a SummarizedExperiment object.
#' Optionally, metadata can be provided, which will be loaded into the `colData` of the SummarizedExperiment object.
#'
#' @param file_path Character. Path to the Sylph query or profile output file (tab-separated format).
#' @param meta_data data.frame. A tibble or data frame containing sample metadata.
#' The **first column must contain sample names** that match exactly with the sample names in the Sylph output.
#' @param variable Character. Name of the input variable to import, either `Adjusted_ANI` or  `Naive_ANI`. Defaults to `Adjusted_ANI`.
#' @param clean_names Logical. If `TRUE`, file paths will be stripped of their directory path and file extension,
#' leaving only the base file name. Defaults to `TRUE`.
#'
#' @return A SummarizedExperiment object with the following components:
#'   \item{assays}{A matrix containing numeric features such as `Adjusted_ANI`, `Taxonomic_abundance`, `Median_cov`, etc.}
#'   \item{rowData}{DataFrame containing metadata for each sequence, such as `Contig_name` and `Genome_file`.}
#'   \item{colData}{DataFrame containing information for each sample (derived from the `Sample_file` field).}
#'
#' @importFrom data.table data.table
#' @importFrom data.table :=
#' @importFrom data.table .GRP
#' @importFrom Matrix sparseMatrix
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(strainspy)
#' library(SummarizedExperiment)
# Read a Sylph file (query or profile) into a SummarizedExperiment object
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path)
#' # View the SummarizedExperiment
#' se
#' # View the assays (numerical matrix)
#' assay(se)
#' # View the rowData (metadata for contigs)
#' rowData(se)
#' # View the colData (metadata for samples)
#' colData(se)
#' }
#' \dontrun{
#'   # Read a Sylph file (query or profile) with associated metadata into a SummarizedExperiment object
#'   example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#'   example_meta <- readr::read_csv(example_meta_path)
#'   example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#'   se <- read_sylph(example_path, example_meta)
#' }
read_sylph <- function(file_path, meta_data=NULL, variable = "Adjusted_ANI",
                       clean_names = TRUE) {
  
  # Check input argument
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("`file_path` must be a character string specifying the path to the Sylph output file.")
  }
  
  if (!file.exists(file_path)) {
    stop(paste0("The file '", file_path, "' does not exist. Please provide a valid file path."))
  }
  
  # Load a small portion of the file to inspect column names
  preview_col_names <- colnames(data.table::fread(
    file_path,
    nrows = 0
  ))
  
  # Identify if the file is a query or profile output based on the presence of specific columns
  is_profile_output <- all(c(
    "Taxonomic_abundance",
    "Sequence_abundance",
    "True_cov",
    "kmers_reassigned"
  ) %in% preview_col_names)
  
  if (is_profile_output) {
    message("Detected Sylph profile output file.")
    valid_var <- c('Taxonomic_abundance','Sequence_abundance','Adjusted_ANI',
                   'True_cov','Median_cov','Mean_cov_geq1','Naive_ANI')
  } else {
    message("Detected Sylph query output file.")
    valid_var <- c('Adjusted_ANI','Naive_ANI')
  }
  
  if (!variable %in% valid_var) {
    stop(paste0("The variable column '", variable, "' is not valid! Please choose from '",
                paste(valid_var, collapse = "', '"), "'"))
  }
  
  required_columns <- c("Sample_file", "Genome_file", "Contig_name", variable)
  
  # Validate that the data contains the expected columns
  missing_columns <- setdiff(required_columns, preview_col_names)
  if (length(missing_columns) > 0) {
    stop(paste0("The following required columns are missing from the Sylph file: ",
                paste(missing_columns, collapse = ", ")))
  }
  
  # Read the Sylph output file
  sylph_data <- data.table::fread(
    file_path,
    na.strings = c("", "NA"),
    select = required_columns
  )
  
  # Calculate row indices using
  sylph_data[, row_indices := .GRP, by = .(Contig_name, Genome_file)]
  
  # Calculate column indices
  sylph_data[, col_indices := .GRP, by = Sample_file]
  
  # Determine dimensions for the sparse matrix
  n_rows <- max(sylph_data$row_indices)
  n_cols <- max(sylph_data$col_indices)
  
  # Merge metadata if provided
  if (is.null(meta_data)){
    # Generate colData
    col_data <- S4Vectors::DataFrame(unique(sylph_data[, .(Sample_file, col_indices)][order(col_indices)]))
    rownames(col_data) <- col_data[,1]
    col_data[['col_indices']] <- NULL
  } else {
    col_data <- unique(sylph_data[, .(Sample_file, col_indices)][order(col_indices)])
    
    # Extract unique sample names
    meta_samples <- unique(meta_data[[1]])
    
    # Warn about mismatched samples
    if (length(missing_from_meta <- setdiff(col_data$Sample_file, meta_samples)) > 0) {
      stop("The following samples from 'sylph_data' are not in 'meta_data': ", paste(missing_from_meta, collapse = ", "))
    }
    
    if (length(missing_from_sylph <- setdiff(meta_samples, col_data$Sample_file)) > 0) {
      stop("The following samples from 'meta_data' are not in 'sylph_data': ", paste(missing_from_sylph, collapse = ", "))
    }
    
    col_data <- S4Vectors::DataFrame(base::merge(col_data, meta_data,
                                                 by.x = "Sample_file",
                                                 by.y = names(meta_data)[1],
                                                 all.x = TRUE)) # Keeps all Sample_file entries from sylph_data
    
    
    rownames(col_data) <- col_data[["Sample_file"]]
    col_data[['col_indices']] <- NULL
  }
  
  # Extract row metadata (rowData)
  # Generate row_data using unique combinations and indices
  row_data <- sylph_data[!duplicated(row_indices),
                         .(row_indices, Contig_name, Genome_file)][order(row_indices)]
  row_data <- S4Vectors::DataFrame(row_data)
  rownames(row_data) <- row_data$Contig_name
  row_data$row_indices <- NULL
  
  # Deal with file names and extensions
  if (clean_names){
    row_data$Genome_file <- tools::file_path_sans_ext(basename(row_data$Genome_file),
                                                      compression = TRUE)
    row_data$Genome_file <- gsub("^.*(GC[A-Z]_\\d+\\.\\d+).*", "\\1", row_data$Genome_file)
    
    col_data$Sample_file <- tools::file_path_sans_ext(basename(col_data$Sample_file),
                                                      compression = TRUE)
    rownames(col_data) <- tools::file_path_sans_ext(basename(rownames(col_data)),
                                                    compression = TRUE)
  }
  
  # Return the SummarizedExperiment object
  return(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(Matrix::sparseMatrix(
        i = sylph_data[['row_indices']],
        j = sylph_data[['col_indices']],
        x = sylph_data[[variable]],
        dims = c(n_rows, n_cols),
        repr = "R" # Specify row-compressed format
      )),
      rowData = row_data,
      colData = col_data
    )
  )
}
