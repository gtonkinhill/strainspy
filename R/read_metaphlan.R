#' read_metaphlan
#' 
#' Read a merged `MetaPhlAn` profile output and create a SummarizedExperiment object. 
#' Optionally, metadata can be provided, which will be loaded into the `colData` of the SummarizedExperiment object.
#'
#' @note Since `MetaPhlAn` profiles are usually generated per sample, first use `merge_metaphlan_files()` to merge and generate the input for this function.
#' 
#' @param file_path Character. Path to the merged metaphlan output file.
#' @param meta_data data.frame. A tibble or data frame containing sample metadata. Metadata requirements are identical to `read_sylph()`.
#' @param variable Character. Name of the input variable to import. Defaults to `relative_abundance`.
#' @param clean_names Logical. If `TRUE`, file paths will be stripped of their directory path and file extension,
#' leaving only the base file name. Defaults to `TRUE`.
#'
#' @return A SummarizedExperiment object with the following components:
#'   \item{assays}{A matrix containing numeric features such as `relative_abundance`}
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
#'   # Read a merged metaphlan profile file into a SummarizedExperiment object
#'   example_path <- system.file("extdata", "metaphlan_merged.tsv.gz", package = "strainspy")
#'   mp <- read_metaphlan(example_path)
#'   # View the SummarizedExperiment
#'   mp
#'   # View the assays (numerical matrix)
#'   assay(mp)
#'   # View the rowData (metadata for contigs)
#'   rowData(mp)
#'   # View the colData (metadata for samples)
#'   colData(mp)
#' }
#' @export
read_metaphlan <- function(file_path, meta_data=NULL, variable="relative_abundance",
                           clean_names = TRUE) {
  
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("`file_path` must be a character string specifying the path to the metaphlan merged output file.")
  }
  
  if (!file.exists(file_path)) {
    stop(paste0("The file '", file_path, "' does not exist. Please provide a valid file path."))
  }
  
  # Load a small portion of the file to inspect column names
  preview_col_names <- colnames(data.table::fread(
    file_path,
    nrows = 0
  ))
  
  required_columns <- c("query_name", "name", variable)
  genome_and_contig = F
  
  
  # Validate that the data contains the expected columns 
  missing_columns <- setdiff(required_columns, preview_col_names)
  if (length(missing_columns) > 0) {
    stop(paste0("The following required columns are missing from the metaphlan file: ",
                paste(missing_columns, collapse = ", ")))
  }
  
  
  # Read the metaphlan output file
  metaphlan_data <- data.table::fread(
    file_path,
    na.strings = c("", "NA"),
    header = T,
    select = required_columns
  )
  
  colnames(metaphlan_data)[1] = "Sample_file" # we'll rename this to match with sylph
  colnames(metaphlan_data)[2] = "Contig_name" # we'll rename this to match with sylph
  
  # Calculate row indices using
  metaphlan_data[, row_indices := .GRP, by = Contig_name]
  
  # Calculate column indices
  metaphlan_data[, col_indices := .GRP, by = Sample_file]
  
  # Determine dimensions for the sparse matrix
  n_rows <- max(metaphlan_data$row_indices)
  n_cols <- max(metaphlan_data$col_indices)
  
  # Merge metadata if provided
  if (is.null(meta_data)){
    # Generate colData
    col_data <- S4Vectors::DataFrame(unique(metaphlan_data[, .(Sample_file, col_indices)][order(col_indices)]))
    rownames(col_data) <- col_data[,1]
    col_data[['col_indices']] <- NULL
  } else {
    col_data <- unique(metaphlan_data[, .(Sample_file, col_indices)][order(col_indices)])
    
    # Extract unique sample names
    meta_samples <- unique(meta_data[[1]])
    
    # Warn about mismatched samples
    if (length(missing_from_meta <- setdiff(col_data$Sample_file, meta_samples)) > 0) {
      stop("The following samples from 'metaphlan_data' are not in 'meta_data': ", paste(missing_from_meta, collapse = ", "))
    }
    
    if (length(missing_from_metaphlan <- setdiff(meta_samples, col_data$Sample_file)) > 0) {
      stop("The following samples from 'meta_data' are not in 'metaphlan_data': ", paste(missing_from_metaphlan, collapse = ", "))
    }
    
    col_data <- S4Vectors::DataFrame(base::merge(col_data, meta_data,
                                                 by.x = "Sample_file",
                                                 by.y = names(meta_data)[1],
                                                 all.x = TRUE)) # Keeps all Sample_file entries from metaphlan_data
    
    
    rownames(col_data) <- col_data[["Sample_file"]]
    col_data[['col_indices']] <- NULL
  }
  
  # Extract row metadata (rowData)
  # Generate row_data using unique combinations and indices
  
  row_data <- metaphlan_data[!duplicated(row_indices),
                             .(row_indices, Contig_name)][order(row_indices)]
  
  
  row_data <- S4Vectors::DataFrame(row_data)
  rownames(row_data) <- row_data$Contig_name
  row_data$row_indices <- NULL
  
  # Deal with file names and extensions
  if (clean_names){
    col_data$Sample_file <- tools::file_path_sans_ext(basename(col_data$Sample_file),
                                                      compression = TRUE)
    rownames(col_data) <- tools::file_path_sans_ext(basename(rownames(col_data)),
                                                    compression = TRUE)
  }
  
  # Return the SummarizedExperiment object
  return(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(Matrix::sparseMatrix(
        i = metaphlan_data[['row_indices']],
        j = metaphlan_data[['col_indices']],
        x = metaphlan_data[[variable]],
        dims = c(n_rows, n_cols),
        repr = "R" # Specify row-compressed format
      )),
      rowData = row_data,
      colData = col_data
    )
  )
  
}

#' merge_metaphlan_files
#' 
#' Merge `MetaPhlAn >=4.0` profiles into a single file to use with `read_metaphlan()` and generate the taxonomy file required for downstream analysis.
#' 
#' @param metaphlan_files A character vector of file paths to merge. Optionally, if the input is named, `names(metaphlan_files)` will be used as sample identifiers in place of the filename.
#' @param output_folder Folder path to save outputs.
#' @param output_filename File name for the two outputs. By default, the two files will be saved as metaphlan_merged.tsv.gz and metaphlan_taxonomy.tsv.gz
#' 
#' @export
merge_metaphlan_files <- function(metaphlan_files, output_folder, output_filename = "metaphlan"){
  if (length(metaphlan_files) == 0) {
    stop("No metaphlan profiles provided.")
  }
  
  header <- data.table::fread(metaphlan_files[1], nrows = 10, header = T)  ### There are some leading comment lines, best to read beyond them nrows = 10
  header = colnames(header)
  
  if(!"#clade_name" %in% header  | !"relative_abundance" %in% header){
    stop('#clade_name and relative_abundance columns are required!')
  }
  
  if(is.null(names(metaphlan_files))){
    xx <- tools::file_path_sans_ext(basename(metaphlan_files),
                                    compression = TRUE)
    xx <- gsub("^.*(GC[A-Z]_\\d+\\.\\d+).*", "\\1", xx)
    names(metaphlan_files) = xx
  }
  
  data_list <- lapply(seq_along(metaphlan_files), function(i) {
    file = metaphlan_files[i]
    tmp = data.table::fread(file);
    tmp = tmp[, .(`#clade_name`, relative_abundance)]
    tmp = tmp[grepl("t__", `#clade_name`)]
    tmp[, query_name := names(metaphlan_files)[i]]
    
    tmp[, genome := sub(".*t__", "", `#clade_name`)]  # Extract everything after "t__"
    tmp[, taxonomy := sub("\\|t__.*", "", `#clade_name`)]  # Remove everything after "|t__"
    tmp[, taxonomy := sub("^k__", "d__", taxonomy)]
    tmp[, taxonomy := gsub("\\|", ";", taxonomy)]
    
  }
  )
  merged_data <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
  
  output_table <- merged_data[, .(query_name, name = genome, relative_abundance = relative_abundance)]
  taxonomy_table <- unique(merged_data[, .(genome, taxonomy)])
  
  data.table::fwrite(taxonomy_table, file = file.path(output_folder, paste0(output_filename, "_taxonomy.tsv.gz")), sep = "\t")
  data.table::fwrite(output_table, file = file.path(output_folder, paste0(output_filename, "_merged.tsv.gz")), sep = "\t")
  
  message("Outputs saved to: ", output_folder, "... You can use these with read_metaphlan()")
}
