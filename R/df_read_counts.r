#' df_read_counts
#'
#' Function to read allelic counts
#'
#' @param file_path path to the file
#' @return df contains allele counts 
#' @examples
#' df_read_counts(file_path) 

df_read_counts <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  read_tsv(file_path, col_types = cols(
    gene_id = col_character(),
    a_count = col_double(),
    b_count = col_double()
  ))
}