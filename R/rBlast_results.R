#' Local Blast Results
#'
#' @param results_table A dataframe of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with Rsamtools::bamFile
#' @param num_results A number indicating number of Metascope results to blast
#' @param num_reads_per_result A number indicating number of reads to blast per result
#' @param hit_list A number of how many blast results to fetch for each read
#' @param db_path Blast database path
#' @param out_path Output directory to save csv files, including base name of files
#'
#' @return Creates and exports num_results number of csv files with blast results from local blast

rBlast_results <- function(results_table, bam_file, num_results = 10, num_reads_per_result = 100, hit_list = "10",
                           db_path, out_path, sample_name = NULL) {
  for (i in seq.int(num_results)) {
    df <- rBLAST_single_result(results_table, bam_file, which_result = i,
                               num_reads = num_reads_per_result, hit_list = hit_list,
                               db_path = db_path)
    tax_id <- result_table[i,1]
    write.csv(df, paste0(out_path, sample_name, "_", "tax_id_", tax_id, "_", i, ".csv"))
  }
}

