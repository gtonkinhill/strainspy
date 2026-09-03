#' Read Sylph query or profile output file and create a SummarizedExperiment object
#'
#' This function reads a Sylph output file (either query or profile), parses the data, and returns a SummarizedExperiment object.
#' An optional metadata table can be provided and added directly to `colData` via an internal call to `modify_metadata()`.
#'
#' If metadata is provided, it must meet the input requirements described in 
#' `modify_metadata()`. If appending metadata fails due to pre processing requirements, 
#' this function will issue a detailed warning and still return the `SummarizedExperiment` 
#' object without metadata. In that case, apply necessary fixes and subsequently 
#' call `modify_metadata()`. 
#'
#' @param file_path Character. Path to the Sylph query or profile output file(s) in `tsv` format. 
#' A single merged sylph output file processed with `merge_sylph_files()`, or a character vector of paths 
#' to individual sylph output files can be provided. If multiple files are provided, they will be merged internally and saved to a temporary file. 
#' If the internal merging fails, use `merge_sylph_files()` to merge the files externally and provide the merged file path.
#' @param meta_data data.frame. An optional tibble or data frame containing sample metadata. See Details
#' @param variable Character. Name of the input variable to import. Defaults to `Adjusted_ANI`, other options can include `Naive_ANI`, `Taxonomic_abundance` and `Sequence_abundance`.
#' @param min_identity Numeric. Minimum identity threshold for filtering ANI values. Defaults to `95`.
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
#' 
#'   # Read a Sylph file 
#'   example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", 
#'   package = "strainspy")
#'   sy <- read_sylph(example_path)
read_sylph <- function(file_path, meta_data=NULL, variable = "Adjusted_ANI", min_identity = 95, clean_names = TRUE) {
  
  # Check input argument
  if(file_path %>% length() > 1){
    message("Paths to multiple Sylph files provided. Merging them into a single file for processing.")
    temp_file <- tempfile(fileext = ".tsv") # We are using disk IO for quicker integration for now, change later
    merge_sylph_files(file_path, temp_file)
    file_path <- temp_file
  }
  
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("`file_path` must be a character string specifying the path to the Sylph output file.")
  }
  
  if (!file.exists(file_path)) {
    stop("The file '", file_path, "' does not exist. Please provide a valid file path.")
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
  
  # Filter ANI values based on the minimum identity threshold
  if (grepl("ANI", variable)) {
    if (!is.numeric(min_identity) || min_identity < 0 || min_identity > 100) {
      stop("`min_identity` must be a numeric value between 0 and 100.")
    }
    sylph_data <- sylph_data[get(variable) >= min_identity]
  }
  
  # Calculate row indices using
  sylph_data[, row_indices := .GRP, by = .(Contig_name, Genome_file)]
  
  # Calculate column indices
  sylph_data[, col_indices := .GRP, by = Sample_file]
  
  # Determine dimensions for the sparse matrix
  n_rows <- max(sylph_data$row_indices)
  n_cols <- max(sylph_data$col_indices)
  
  # Merge metadata if provided
  # if (is.null(meta_data)){
  # Generate colData
  col_data <- S4Vectors::DataFrame(unique(sylph_data[, .(Sample_file, col_indices)][order(col_indices)]))
  rownames(col_data) <- col_data[,1]
  col_data[['col_indices']] <- NULL
  # } else {
  # col_data <- unique(sylph_data[, .(Sample_file, col_indices)][order(col_indices)])
  # 
  # # Extract unique sample names
  # meta_samples <- unique(meta_data[[1]])
  # 
  # # Warn about mismatched samples
  # if (length(missing_from_meta <- setdiff(col_data$Sample_file, meta_samples)) > 0) {
  #   stop("The following samples from 'sylph_data' are not in 'meta_data': ", paste(missing_from_meta, collapse = ", "))
  # }
  # 
  # if (length(missing_from_sylph <- setdiff(meta_samples, col_data$Sample_file)) > 0) {
  #   stop("The following samples from 'meta_data' are not in 'sylph_data': ", paste(missing_from_sylph, collapse = ", "))
  # }
  # 
  # col_data <- S4Vectors::DataFrame(base::merge(col_data, meta_data,
  #                                              by.x = "Sample_file",
  #                                              by.y = names(meta_data)[1],
  #                                              all.x = TRUE)) # Keeps all Sample_file entries from sylph_data
  # 
  # 
  # rownames(col_data) <- col_data[["Sample_file"]]
  # col_data[['col_indices']] <- NULL
  # }
  
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
  
  # Duplicated names would silently misalign metadata downstream.
  check_unique_names(rownames(col_data), "sample", clean_names)
  check_unique_names(rownames(row_data), "feature", clean_names)

  se = SummarizedExperiment::SummarizedExperiment(
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
  
  if(!is.null(meta_data)){ # User has provided meta_data, attempt to automatically append it to se
    se <- tryCatch({
      modify_metadata(se, meta_data, replace = TRUE)
    }, error = function(e) {
      warning("Automated attachment of metadata failed: ", conditionMessage(e), 
              "\nReturning SummarizedExperiment without metadata. Use `modify_metadata()` after applying necessary fixes.")
      se  # return original se
    })
  }
  
  # Return the SummarizedExperiment object
  return(se)
}

#' merge_sylph_files
#' 
#' Merge the outputs of `sylph query` or `sylph profile` into a single file.
#' The output of this function is suitable to use in read_sylph()
#' 
#' 
#' @param sylph_files A character vector of sylph output file paths to merge.
#' @param output Save path for the output file. Internally, data.table::fwrite is used, provide the extension `gz` to save a compressed file.
#' @return Invisibly returns `NULL`. The merged sylph table is written to `output`.
#' @examples
#' # Split the bundled profile into two per-sample files, then merge them back.
#' profile <- data.table::fread(
#'   system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' )
#' samples <- head(unique(profile$Sample_file), 2)
#' parts <- vapply(samples, function(s) {
#'   f <- tempfile(fileext = ".tsv")
#'   data.table::fwrite(profile[profile$Sample_file == s, ], f, sep = "\t")
#'   f
#' }, character(1))
#'
#' merged <- tempfile(fileext = ".tsv.gz")
#' merge_sylph_files(parts, merged)
#'
#' # The merged file carries one header and both samples
#' se <- read_sylph(merged)
#' dim(se)
#' @export
merge_sylph_files <- function(sylph_files, output) {
  stopifnot(is.character(sylph_files), length(sylph_files) > 0)
  stopifnot(is.character(output), length(output) == 1)
  
  if (!all(file.exists(sylph_files))) {
    missing <- sylph_files[!file.exists(sylph_files)]
    stop("Sylph files not found: ", paste(missing, collapse = ", "))
  }
  
  open_input <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, open = "rt")
    } else {
      file(path, open = "rt")
    }
  }
  
  open_output <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, open = "wt")
    } else {
      file(path, open = "wt")
    }
  }
  
  out <- open_output(output)
  on.exit(close(out), add = TRUE)
  
  for (i in seq_along(sylph_files)) {
    con <- open_input(sylph_files[i])
    
    if (i > 1) {
      readLines(con, n = 1)  # skip header
    }
    
    repeat {
      lines <- readLines(con, n = 10000)
      if (!length(lines)) break
      writeLines(lines, out)
    }
    
    close(con)
  }
  
  invisible(output)
}
