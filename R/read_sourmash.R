#' read_sourmash
#' 
#' This function reads a `sourmash gather` or `sourmash search` -like output CSV file and create a SummarizedExperiment object. 
#' An optional metadata table can be provided and added directly to `colData` via an internal call to `modify_metadata()`.
#'
#' @note Some of sourmash subcommands like `gather` are designed to take exactly one query file. 
#' To merge and save a CSV file suitable to use with this function, use `merge_sourmash_files()`
#' #' If metadata is provided, it must meet the input requirements described in 
#' `modify_metadata()`. If appending metadata fails due to pre processing requirements, 
#' this function will issue a detailed warning and still return the `SummarizedExperiment` 
#' object without metadata. In that case, apply necessary fixes and subsequently 
#' call `modify_metadata()`. 
#' 
#' @param file_path Character. Path to the sourmash gather file (CSV format).
#' @param meta_data data.frame. A tibble or data frame containing sample metadata. Metadata requirements are identical to `read_sylph()`.
#' @param variable Character. Name of the input variable to import. Defaults to `match_containment_ani`.
#' @param clean_names Logical. If `TRUE`, file paths will be stripped of their directory path and file extension,
#' leaving only the base file name. Defaults to `TRUE`.
#'
#' @return A SummarizedExperiment object with the following components:
#'   \item{assays}{A matrix containing numeric features such as `Adjusted_ANI`, `Taxonomic_abundance`, `Median_cov`, etc.}
#'   \item{rowData}{DataFrame containing metadata for each sequence, such as `Contig_name` and `Genome_file`.}
#'   \item{colData}{DataFrame containing information for each sample (derived from the `Sample_file` field).}
#'
#' @importFrom data.table data.table merge.data.table
#' @importFrom data.table :=
#' @importFrom data.table .GRP
#' @importFrom Matrix sparseMatrix
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame na.omit
#'
#' @examples
#' \dontrun{
#'   # Read a merged sourmash file into a SummarizedExperiment object
#'   example_path <- system.file("extdata", "example_sourmash_merged.csv.gz", package = "strainspy")
#'   sm <- read_sourmash(example_path)
#'   # View the SummarizedExperiment
#'   sm
#'   # View the assays (numerical matrix)
#'   assay(sm)
#'   # View the rowData (metadata for contigs)
#'   rowData(sm)
#'   # View the colData (metadata for samples)
#'   colData(sm)
#' }
#' @export
read_sourmash <- function(file_path, meta_data=NULL, variable="match_containment_ani",
                          clean_names = TRUE) {
  
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("`file_path` must be a character string specifying the path to the Sourmash output file.")
  }
  
  if (!file.exists(file_path)) {
    stop(paste0("The file '", file_path, "' does not exist. Please provide a valid file path."))
  }
  
  # Load a small portion of the file to inspect column names
  preview_col_names <- colnames(data.table::fread(
    file_path,
    nrows = 0
  ))
  
  # name is somewhat arbitrary in sourmash. If processed using merge_sourmash_files(),
  # both name = genome_name and contig will be present
  
  if("contig" %in% preview_col_names) {
    required_columns = c("query_name", "name", "contig", variable)
    genome_and_contig = T
  } else {
    required_columns <- c("query_name", "name", variable)
    genome_and_contig = F
  }
  
  # Validate that the data contains the expected columns
  missing_columns <- setdiff(required_columns, preview_col_names)
  if (length(missing_columns) > 0) {
    stop(paste0("The following required columns are missing from the Sourmash file: ",
                paste(missing_columns, collapse = ", ")))
  }
  
  
  # Read the Sourmash output file
  sourmash_data <- data.table::fread(
    file_path,
    na.strings = c("", "NA"),
    header = T,
    select = required_columns
  )
  
  if(genome_and_contig){
    colnames(sourmash_data)[1] = "Sample_file" # we'll rename this to match with sylph
    colnames(sourmash_data)[2] = "Genome_file" # we'll rename this to match with sylph
    colnames(sourmash_data)[3] = "Contig_name" # we'll rename this to match with sylph
  } else {
    colnames(sourmash_data)[1] = "Sample_file" # we'll rename this to match with sylph
    colnames(sourmash_data)[2] = "Contig_name" # we'll rename this to match with sylph
  }
  
  # Calculate row indices using
  sourmash_data[, row_indices := .GRP, by = Contig_name]
  
  # Calculate column indices
  sourmash_data[, col_indices := .GRP, by = Sample_file]
  
  # Determine dimensions for the sparse matrix
  n_rows <- max(sourmash_data$row_indices)
  n_cols <- max(sourmash_data$col_indices)
  
  # Merge metadata if provided
  # if (is.null(meta_data)){
  # Generate colData
  col_data <- S4Vectors::DataFrame(unique(sourmash_data[, .(Sample_file, col_indices)][order(col_indices)]))
  rownames(col_data) <- col_data[,1]
  col_data[['col_indices']] <- NULL
  # } else {
  #   col_data <- unique(sourmash_data[, .(Sample_file, col_indices)][order(col_indices)])
  # 
  #   # Extract unique sample names
  #   meta_samples <- unique(meta_data[[1]])
  # 
  #   # Warn about mismatched samples
  #   if (length(missing_from_meta <- setdiff(col_data$Sample_file, meta_samples)) > 0) {
  #     stop("The following samples from 'sourmash_data' are not in 'meta_data': ", paste(missing_from_meta, collapse = ", "))
  #   }
  # 
  #   if (length(missing_from_sourmash <- setdiff(meta_samples, col_data$Sample_file)) > 0) {
  #     stop("The following samples from 'meta_data' are not in 'sourmash_data': ", paste(missing_from_sourmash, collapse = ", "))
  #   }
  # 
  #   col_data <- S4Vectors::DataFrame(base::merge(col_data, meta_data,
  #                              by.x = "Sample_file",
  #                              by.y = names(meta_data)[1],
  #                              all.x = TRUE)) # Keeps all Sample_file entries from sourmash_data
  # 
  # 
  #   rownames(col_data) <- col_data[["Sample_file"]]
  #   col_data[['col_indices']] <- NULL
  # }
  
  # Extract row metadata (rowData)
  # Generate row_data using unique combinations and indices
  
  if(genome_and_contig){
    row_data <- sourmash_data[!duplicated(row_indices),
                              .(row_indices, Contig_name, Genome_file)][order(row_indices)]
  } else {
    row_data <- sourmash_data[!duplicated(row_indices),
                              .(row_indices, Contig_name)][order(row_indices)]
  }
  
  row_data <- S4Vectors::DataFrame(row_data)
  rownames(row_data) <- row_data$Contig_name
  row_data$row_indices <- NULL
  
  # Deal with file names and extensions
  if (clean_names){
    if(genome_and_contig){
      row_data$Genome_file <- tools::file_path_sans_ext(basename(row_data$Genome_file),
                                                        compression = TRUE)
      row_data$Genome_file <- gsub("^.*(GC[A-Z]_\\d+\\.\\d+).*", "\\1", row_data$Genome_file)
    }
    
    col_data$Sample_file <- tools::file_path_sans_ext(basename(col_data$Sample_file),
                                                      compression = TRUE)
    rownames(col_data) <- tools::file_path_sans_ext(basename(rownames(col_data)),
                                                    compression = TRUE)
  }
  
  se = SummarizedExperiment::SummarizedExperiment(
    assays = list(Matrix::sparseMatrix(
      i = sourmash_data[['row_indices']],
      j = sourmash_data[['col_indices']],
      x = sourmash_data[[variable]],
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

#' merge_sourmash_files
#' 
#' Merge the outputs of `sourmash search` or `sourmash gather` into a single file.
#' The output of this function is suitable to use in read_sourmash()
#' 
#' 
#' Refer to \href{https://sourmash.readthedocs.io/en/latest/classifying-signatures.html#searching-for-similar-samples-with-search}{sourmash documentation} for more information.
#' 
#' 
#' @param sourmash_files A character vector of file paths to merge.
#' @param manifest_file If one or several indexed databases were used for searching or gathering, required genome column can be missing.
#' The path to the output of `sourmash sig manifest` can be provided to try and add back the data automatically (Default = NULL). 
#' If provided, columns `md5`, `name` and `filename` are required.  This function uses `name` and `md5` to match the manifest data with sourmash outputs. 
#' In the final output, `filename` will be used to identify each hit. If the `filename` is a path, the basename without extensions will be used instead.
#' @param output Save path for the output file. Internally, data.table::fwrite is used, provide the extension `gz` to save a compressed file.
#' @param strip_unusable Remove columns containing data not useful for further analysis. Default T. 
#' @export
merge_sourmash_files <- function(sourmash_files, manifest_file = NULL, output, strip_unusable = T){
  if (length(sourmash_files) == 0) {
    stop("No sourmash output files provided.")
  }
  
  header <- data.table::fread(sourmash_files[1], nrows = 1, header = T)
  header = colnames(header)
  
  data_list <- lapply(sourmash_files, function(file) fread(file, skip = 1, quote = "\"" ))
  merged_data <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
  
  # if filename is in header, let's name it filename_db
  if("filename" %in% header) header[header == "filename"] = "filename_db"
  
  colnames(merged_data) = header
  
  if(!is.null(manifest_file)){ # we only want the genome name for taxonomic stuff...
    manifest = data.table::fread(manifest_file)
    required_columns <- c("md5", "name", "filename")
    missing_columns <- setdiff(required_columns, colnames(manifest))
    if(length(missing_columns) == 0){
      manifest = data.table(md5 = manifest$md5, name = manifest$name, filename = manifest$filename)
    } else {
      stop("Manifest file doesn't contain all required columns")
    }
    # content in md5 and name must match with md5 and name in merged_data
    if(manifest$md5[1] %in% merged_data$md5 & manifest$name[1] %in% merged_data$name) {
      merged_data = merge.data.table(merged_data, manifest, by =  c("md5", "name"))
    }
    
    # change filename to name and name to contig
    colnames(merged_data)[which("name" == colnames(merged_data))] = "contig"
    colnames(merged_data)[which("filename" == colnames(merged_data))] = "name"
    
    merged_data$name = unname(sapply(basename(merged_data$name),  function(x) sub("\\.[^.]+$", "", sub("\\.[^.]+$", "", x))))
    
  }
  
  if(strip_unusable) {
    drp = c()
    # md5s 
    drp_tmp = grep("md5", colnames(merged_data))
    if(length(drp_tmp) > 0) {
      cat("Stripping", length(drp_tmp), "md5 columns:", paste(colnames(merged_data)[drp_tmp], collapse = ', '), "\n")
      drp = c(drp_tmp)
    }
    
    # NAs
    drp_tmp <- which(sapply(merged_data, function(x) all(is.na(x))))
    if(length(drp_tmp) > 0) {
      cat("Stripping", length(drp_tmp), "all NA columns:", paste(colnames(merged_data)[drp_tmp], collapse = ', '), "\n")
      drp = c(drp, drp_tmp)
    }
    
    # Invariant columns
    drp_tmp <- which(sapply(merged_data, function(x) length(na.omit(unique(x))) == 1))
    if(length(drp_tmp) > 0) {
      cat("Stripping", length(drp_tmp), "invariant columns:", paste(colnames(merged_data)[drp_tmp], collapse = ', '), "\n")
      drp = c(drp, drp_tmp)
    }
    
    set(merged_data, j = colnames(merged_data)[unique(drp)], value = NULL)
  }
  
  data.table::fwrite(merged_data, output, quote = TRUE, col.names = T)
  
  message("Merged sourmash output file saved to: ", output, "... You can use this output with read_sourmash()")
}
