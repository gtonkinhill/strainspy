#' Extract data from a StrainSpy `SummarizedExperiment` object in long format
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#'
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param variables One or multiple variables in `colnames(se@colData)`, except `Sample_File`.
#' @param contigs A vector of contigs to generate violin plots
#' @param drop_zeros Bool. If TRUE, ANI or Abundance 0 values will not be included in the violin plot. Default FALSE.
#'
#' @examples
#' \dontrun{
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
#' example_meta <- readr::read_csv(example_meta_path)
#' # Change variable from character to factor
#' example_meta$Case_status = factor(example_meta$Case_status)
#'
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#'
#' tax_path <- system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy")
#' tax <- read_taxonomy(tax_path)
#'
#' extract_strainspy(se, variables = c("Case_status", "Age_at_collection"), contigs = rownames(se)[1:5], taxonomy=tax)
#'
#' }
#'
#' @export
extract_strainspy <- function(se, variables=NULL, contigs=NULL, taxonomy=NULL, drop_zeros = FALSE){

  # default variables to colnames of colData
  if (is.null(variables)) {
    variables = colnames(SummarizedExperiment::colData(se))
  }

  if(! all(variables %in% colnames(SummarizedExperiment::colData(se))) ){
    stop("At least one variable not found in colnames(se@colData)")
  }

  # check contigs
  chk_contigs = contigs %in% rownames(se)
  if( ! all(chk_contigs)  ){
    cat(paste(contigs[!chk_contigs], '\n', sep = ""))
    stop("Above contigs are missing from rownames(se)")
  }

  tmp <- tibble::rownames_to_column(
    as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ])))),
    var = "sample")
  tmp <- dplyr::bind_cols(as_tibble(se@colData[variables]), tmp)

  df_long <- tidyr::pivot_longer(tmp, cols = -c(sample, dplyr::all_of(variables)), names_to = "Contig_name", values_to = "ANI")

  if (drop_zeros) {
    df_long <- df_long |> filter(ANI > 0)
  }

  # Add other columns from rowdata
  df_long <- df_long |>
    left_join(tibble::as_tibble(SummarizedExperiment::rowData(se)), by = "Contig_name")

  # TODO: add some checks and informative error messages
  if (!is.null(taxonomy)) {
    df_long <- df_long |>
      left_join(taxonomy, by = c("Genome_file" = "Genome"))
  }

  # Reorder columns in df_long
  df_long <- df_long |>
    dplyr::select(sample, dplyr::all_of(variables), Contig_name, dplyr::everything(), ANI)

  return(df_long)
}
