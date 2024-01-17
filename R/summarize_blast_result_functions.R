#'
#'

summarize_blast_results <- function(result_table, blast_results_dir){
  all_files <- list.files(path = blast_results_dir, full.names = TRUE)
  final_blast_result <- list()
  for (current_path in all_files){
    tryCatch(
      {
        current_blast_result <- read.csv(current_path, header = TRUE)
        blast_result_summary <- current_blast_result %>%
          group_by(qseqid) %>%
          summarise(n = n(), evalue = mean(evalue)) %>%
          top_n(1,n)
        final_blast_result <- rbind(final_blast_result, blast_result_summary)
        },
    error = function(e) {
      message(conditionMessage(e))
    })
  }
  return(final_blast_result)
}

test_blast_result <- read.csv(all_files[1],header = TRUE)
test_blast_summary <- list()
for (query in unique(test_blast_result$qseqid)){
  current_test_blast_result <- test_blast_result[test_blast_result$qseqid == query, ]
  total_hits <- nrow(current_test_blast_result)
  current_test_blast_summary <- current_test_blast_result %>%
    group_by(name) %>%
    summarise(queryid = query, n = n(), frac = n()/total_hits, evalue = mean(evalue))
  test_blast_summary <- rbind(test_blast_summary, current_test_blast_summary[!(current_test_blast_summary$frac < 0.5),])
}

all_reads_summary <- test_blast_summary %>%
  group_by(name) %>%
  summarise(n = n())

uniqueness_score <- 1/nrow(all_reads_summary)





## testing with bls11_S46
bamFile <- Rsamtools::BamFile("E:/bls11_S46.bam")
seq_info_df <- data.frame(Rsamtools::seqinfo(bamFile))
seq_info_df$seqnames <- row.names(seq_info_df)
allGenomes <- grep("ti|652616", seq_info_df$seqnames, value = TRUE, fixed = TRUE)
results_df <- read.csv("E:/bls11_S46.metascope_id.csv", header= TRUE)


# Getting blast metrics and appending to metascope results
results_df <- read.csv("E:/bls1_S36.csv", header = TRUE, skip = 1)
blast_results_df <- read.csv("E:/blast_results/bls1_S36_blast_1.csv", header = TRUE)
blast_result_metrics(blast_results_df)
blast_files <- list.files("E:/blast_results/", full.names = TRUE)
blast_result_summary <- data.frame(best_hit = c(), uniqueness_score = c(), percentage_hit = c(), contaminant_score = c())
for (i in 1:10) {
  blast_file <- paste0("E:/blast_results/bls1_S36_blast_", i, ".csv")
  print(blast_result_metrics(blast_file))
  blast_result_summary <- rbind(blast_result_summary, blast_result_metrics(blast_file))
  cat("current file", blast_file)
}
View(cbind(results_df[1:10,], blast_result_summary))
