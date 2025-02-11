#' process_samples
#'
#' Main function to process multiple samples
#'
#' @param samples Description of the first parameter
#' @param input_rna input RNA df
#' @param input_dna input DNA df
#' @param output_dir output path
#' @param z
#' @param ABminCut
#' @param ABminCut_noCI
#' @param Threshold
#' @return the output print of IB analysis
#' @examples
#' process_samples(samples, input_rna, input_dna, output_dir, z=1.96, ABminCut=1, ABminCut_noCI=0.6)

process_samples <- function(samples, input_rna, input_dna, output_dir, z=1.96, ABminCut=1, ABminCut_noCI=0.6,
                            Threshold = 10) {

  Z <- z
  ABminCuT <- ABminCut
  ABminCuT_noCI <- ABminCut_noCI
  ThresholD <- Threshold

  for (sample in samples) {
    message("Processing: ", sample)

    file_rna <- file.path(input_rna, paste0(sample, "_RNA_counts.txt"))
    file_dna <- file.path(input_dna, paste0(sample, "_DNA_counts.txt"))

    df_rna <- df_read_counts(file_rna)
    df_dna <- df_read_counts(file_dna)

    results <- calculate_ABratio(df_rna, df_dna, z=Z, ABminCut=ABminCuT, ABminCut_noCI=ABminCuT_noCI, Threshold=ThresholD)

    if (!is.null(results)) {
      write_tsv(results, file.path(output_dir, paste0(sample, "_AI_results.txt")))
    }
  }
}
