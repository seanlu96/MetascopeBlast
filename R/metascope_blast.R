#' MetaScope Blast
#'
#' @param metascope_id_path
#'
#' @export
#'

metascope_blast <- function(metascope_id_path, tmp_dir, out_dir, sample_name,
                            num_reads = 100, hit_list = "10", db_path,
                            uniqueness_score_by = "species",
                            percentage_hit_by = "species",
                            contaminant_score_by = "genus") {
  # Sort and index bam file
  bam_file_path <- list.files(path = tmp_dir, full.names = TRUE)
  sorted_bam_file_path <- file.path(tmp_dir, paste0(sample_name, "_sorted"))
  Rsamtools::sortBam(bam_file_path, destination = sorted_bam_file_path)
  sorted_bam_file <- paste0(sorted_bam_file_path, ".bam")
  Rsamtools::indexBam(sorted_bam_file)
  bam_file <- Rsamtools::BamFile(sorted_bam_file, index = sorted_bam_file)

  # Load in metascope id file and clean unknown genomes
  metascope_id <- read.csv(metascope_id_path, header = TRUE)

  # Create blast directory in tmp directory to save blast results in
  blast_tmp_dir <- dir.create(file.path(tmp_dir,"blast"))

  # Run rBlast on all metascope microbes
  rBlast_results(metascope_id, bam_file, num_results = nrow(metascope_id),
                 num_reads_per_result = num_reads, hit_list = hit_list,
                 db_path = db_path, out_path = blast_tmp_dir, sample_name = sample_name)

  # Run Blast metrics
  blast_result_metrics_list <- lapply(list.files(blast_tmp_dir, full.names = TRUE),
                                      blast_result_metrics(),
                                      uniqueness_score_by = uniqueness_score_by,
                                      percentage_hit_by = percentage_hit_by,
                                      contaminant_score_by = contaminant_score_by)

  # Append Blast Metrics to MetaScope results
  blast_result_metrics_df <- as.data.frame(do.call(rbind, blast_result_metrics_list))
  metascope_blast_df <- cbind(metascope_id, blast_result_metrics_df)
  write.csv(file.path(out_dir, paste0(sample_name, ".metascope_blast.csv")))

  # TODO: Apply Filters for blast metrics
}
