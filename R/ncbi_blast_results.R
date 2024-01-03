#' Blast Results to NCBI API
#'
#' @param results_table A dataframe of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with Rsamtools::bamFile
#' @param num_results A number indicating number of Metascope results to blast
#' @param num_reads_per_result A number indicating number of reads to blast per result
#' @param hit_list A number of how many blast results to fetch for each read
#' @param out_path Output directory to save csv files, including basename of files
#'
#' @return Creates and exports num_results number of csv files with blast results from NCBI

ncbi_blast_results <- function(results_table, bam_file, num_results = 10, num_reads_per_result = 10, hit_list = "10",
                               out_path) {
  for (i in seq.int(num_results)) {
    df <- blast_single_result(results_table, bam_file, which_result = i,
                              num_reads = num_reads_per_result, hit_list = hit_list)
    write.csv(df, paste0(out_path, i, ".csv"))
  }
}
