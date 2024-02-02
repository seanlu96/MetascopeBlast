# 1. takes id and bam file, finds all sequence names associated with (it should just be one), maybe just take the first match!
# Get ID and full seqname
# This function only returns a vector of the sequences
getSeqs <- function(id, bamFile, n = 10) {
  # Get sequence info (Genome Name) from Bam file
  seq_info_df <- data.frame(Rsamtools::seqinfo(bamFile))
  seq_info_df$seqnames <- row.names(seq_info_df)
  allGenomes <- grep(paste0("ti|", id), seq_info_df$seqnames, value = TRUE, fixed = TRUE)
  # Sample one of the Genomes that match
  Genome = sample(allGenomes, 1)
  # Scan Bam file for all sequences that match genome
  param <- ScanBamParam(what = c("rname", "seq"),
                        which = GRanges(Genome, IRanges(1, 1e+07)))
  allseqs <- scanBam(bamFile, param = param)[[1]]
  seqs <- as.character(allseqs$seq)
  seqs <- sample(seqs, n)

  #print(x)
  return(seqs)
}


rBLAST_single_result <- function(results_table, bam_file, which_result = 1, num_reads = 100, hit_list = "10", db_path, quiet = TRUE) {
  res <- tryCatch( #If any errors, should just skip the organism
    {
      genome_name <- results_table[which_result,2]
      if (!quiet) message("Current id: ", genome_name)
      tax_id <- results_table[which_result,1]
      if (!quiet) message("Current ti: ", tax_id)
      fasta_seqs <- Biostrings::DNAStringSet(getSeqs(id = tax_id, bamFile = bam_file, n = num_reads))
      blast_db <- blast(db = db_path, type = "blastn")
      df <- predict(blast_db, fasta_seqs,
                    custom_format ="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids",
                    BLAST_args = "-max_target_seqs 10")
      df$MetaScope_Taxid <- tax_id
      df$MetaScope_Genome <- genome_name
      df
    },
    error = function(e) {
      cat("Error", conditionMessage(e))
      df <- data.frame(qseqid=NA, sseqid=NA, pident=NA, length=NA,
                       mismatch=NA, gapopen=NA, qstart=NA, qend=NA,
                       sstart=NA, send=NA, evalue=NA, bitscore=NA)
      df$ti <- tax_id
      df$genome <- genome_name
      df
    }
  )
  return(df)
}





#' rBlast_results
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
    tax_id <- results_table[i,1]
    write.csv(df, file.path(out_path, paste0(sample_name, "_", "tax_id_", tax_id, "_", i, ".csv")))
  }
}

