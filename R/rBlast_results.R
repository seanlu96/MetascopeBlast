rBlast_results <- function(results_table, bam_file, num_results = 10, num_reads_per_result = 100, hit_list = "10",
                           db_path, out_path) {
  for (i in seq.int(num_results)) {
    df <- rBLAST_single_result(results_table, bam_file, which_result = i,
                               num_reads = num_reads_per_result, hit_list = hit_list,
                               db_path = db_path)
    write.csv(df, paste0(out_path, i, ".csv"))
  }
}

