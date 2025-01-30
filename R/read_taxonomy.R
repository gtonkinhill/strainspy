#' Read and process a taxonomy file in Sylph format
#'
#' This function reads a taxonomy file in Sylph format.
#'
#' @param file_path Path to the taxonomy file.
#' @return A data.table containing parsed taxonomy information.
#' @export
#' @importFrom data.table fread tstrsplit
#' @importFrom dplyr arrange
read_taxonomy <- function(file_path) {
  # Read the taxonomy file
  tax <- data.table::fread(file_path, sep = "\t", col.names = c("Genome", "Taxonomy"))

  # Clean and split the Taxonomy column into individual taxonomic ranks
  tax$Taxonomy <- gsub("[a-z]__", "", tax$Taxonomy)
  tax[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") := data.table::tstrsplit(Taxonomy, ";", fill = TRUE)]

  # Add index and arrange
  tax$index <- 1:nrow(tax)
  tax <- tax |>
    dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, Species, index)
  tax$index <- 1:nrow(tax)

  return(tax)
}

# Example usage
# taxonomy_data <- read_taxonomy("./inst/extdata/sylph_gtdb_rep_taxonomy.tsv.gz")
