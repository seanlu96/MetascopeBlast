library(Rsamtools)
library(taxize)
library(httr)
library(xml2)
library(Biostrings)
library(dplyr)
#library(tidyverse)
library(tidyr)
library(R.utils)
library(rBLAST)


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

# Sends FASTA formated string to NCBI BLAST URL API

blastSeq <- function (x, database = "nt", hitListSize = "20",
                      filter = "L", expect = "10", program = "blastn",
                      attempts = 10, baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi") {
  query <- list(
    CMD = "Put",
    QUERY = as.character(x),
    DATABASE = database,
    HITLIST_SIZE = hitListSize,
    FILTER = filter,
    EXPECT = expect,
    PROGRAM = program
  )

  # Make the initial request
  response <- POST(baseUrl, body = query)
  print(status_code(response))


  # Check if the request was successful
  if (status_code(response) == 200) {
    # Extract RID from the response
    x <- content(response, "text")
    #print(paste("This is the content:", x))

    rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
    rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
    print(paste("This is the RID:", rid))
    print(paste("This is the rtoe:", rtoe))

    # Construct the URL for retrieving results
    url1 <- paste0(baseUrl, "?RID=", rid, "&FORMAT_TYPE=XML&CMD=Get")
    print(paste("This is the URL:", url1))

    Sys.sleep(rtoe)

    # Try to parse the result
    result <- tryParseResult(url1, attempts)

    # qseq <- xml_find_all(result, "//Hsp_qseq") %>% xml_text()
    # hseq <- xml_find_all(result, "//Hsp_hseq") %>% xml_text()

    # res <- list()
    # for (i in seq_len(length(qseq))) {
    #  res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]),
    #                                 rowmask = as(IRanges(), "NormalIRanges"),
    #                                 colmask = as(IRanges(), "NormalIRanges"))
    return(result)
  }
  else {
    stop("Failed to make the initial request.")
  }
}

# Helper function to parse XML content with retries from NCBI API
tryParseResult <- function(url, attempts = 30, sleep_time = 60) {
  for (i in 1:(attempts + 1)) {
    result <- tryCatch({
      xml2::read_xml(url, quiet = TRUE)
    }, error=function(err) NULL)
    if (!is.null(result)) return(result)
    Sys.sleep(sleep_time)
  }
  stop(paste("No results after ", attempts,
             " attempts; please try again later", sep = ""))
}


# Reads XML formatted BLAST results and converts to dataframe

blastXmlDf <- function(result) {
  iterations <- xml_find_all(result, "//Iteration")

  # Initialize empty lists to store data
  iter_nums <- list()
  Hit_nums <- list()
  Hit_ids <- list()
  Hit_defs <- list()
  Hit_accessions <- list()
  Hit_lens <- list()
  Hsp_nums <- list()
  Hsp_bit_scores <- list()
  Hsp_scores <- list()
  Hsp_evalues <- list()
  Hsp_query_froms <- list()
  Hsp_query_tos <- list()
  Hsp_hit_froms <- list()
  Hsp_hit_tos <- list()
  Hsp_query_frames <- list()
  Hsp_hit_frames <- list()
  Hsp_identitys <- list()
  Hsp_positives <- list()
  Hsp_gapss <- list()
  Hsp_align_lens <- list()
  Hsp_qseqs <- list()
  Hsp_hseqs <- list()
  Hsp_midlines <- list()

  # Iterate through Iteration nodes
  for (iteration in iterations) {
    # Extract Iteration_iter-num
    iter_num <- xml_text(xml_find_first(iteration, "./Iteration_iter-num"))

    # Extract Hit nodes
    hits <- xml_find_all(iteration, ".//Hit")

    # Iterate through Hit nodes
    for (hit in hits) {
      # Extract Hit_num and Hit_id
      Hit_num = xml_text(xml_find_first(hit, "./Hit_num"))
      Hit_id = xml_text(xml_find_first(hit, "./Hit_id"))
      Hit_def = xml_text(xml_find_first(hit, "./Hit_def"))
      Hit_accession = xml_text(xml_find_first(hit, "./Hit_accession"))
      Hit_len = xml_text(xml_find_first(hit, "./Hit_len"))
      Hsp_num = xml_text(xml_find_first(hit, ".//Hsp_num"))
      Hsp_bit_score = xml_text(xml_find_first(hit, ".//Hsp_bit-score"))
      Hsp_score = xml_text(xml_find_first(hit, ".//Hsp_score"))
      Hsp_evalue = xml_text(xml_find_first(hit, ".//Hsp_evalue"))
      Hsp_query_from = xml_text(xml_find_first(hit, ".//Hsp_query-from"))
      Hsp_query_to = xml_text(xml_find_first(hit, ".//Hsp_query-to"))
      Hsp_hit_from = xml_text(xml_find_first(hit, ".//Hsp_hit-from"))
      Hsp_hit_to = xml_text(xml_find_first(hit, ".//Hsp_hit-to"))
      Hsp_query_frame = xml_text(xml_find_first(hit, ".//Hsp_query-frame"))
      Hsp_hit_frame = xml_text(xml_find_first(hit, ".//Hsp_hit-frame"))
      Hsp_identity = xml_text(xml_find_first(hit, ".//Hsp_identity"))
      Hsp_positive = xml_text(xml_find_first(hit, ".//Hsp_positive"))
      Hsp_gaps = xml_text(xml_find_first(hit, ".//Hsp_gaps"))
      Hsp_align_len = xml_text(xml_find_first(hit, ".//Hsp_align-len"))
      Hsp_qseq = xml_text(xml_find_first(hit, ".//Hsp_qseq"))
      Hsp_hseq = xml_text(xml_find_first(hit, ".//Hsp_hseq"))
      Hsp_midline = xml_text(xml_find_first(hit, ".//Hsp_midline"))

      # Append data to lists
      iter_nums <- c(iter_nums, rep(iter_num, length(Hit_num)))
      Hit_nums <- c(Hit_nums, Hit_num)
      Hit_ids <- c(Hit_ids, Hit_id)
      Hit_defs <- c(Hit_defs, Hit_def)
      Hit_accessions <- c(Hit_accessions, Hit_accession)
      Hit_lens <- c(Hit_lens, Hit_len)
      Hsp_nums <- c(Hsp_nums, Hsp_num)
      Hsp_bit_scores <- c(Hsp_bit_scores, Hsp_bit_score)
      Hsp_scores <- c(Hsp_scores, Hsp_score)
      Hsp_evalues <- c(Hsp_evalues, Hsp_evalue)
      Hsp_query_froms <- c(Hsp_query_froms, Hsp_query_from)
      Hsp_query_tos <- c(Hsp_query_tos, Hsp_query_to)
      Hsp_hit_froms <- c(Hsp_hit_froms, Hsp_hit_from)
      Hsp_hit_tos <- c(Hsp_hit_tos, Hsp_hit_to)
      Hsp_query_frames <- c(Hsp_query_frames, Hsp_query_frame)
      Hsp_hit_frames <- c(Hsp_hit_frames, Hsp_hit_frame)
      Hsp_identitys <- c(Hsp_identitys, Hsp_identity)
      Hsp_positives <- c(Hsp_positives, Hsp_positive)
      Hsp_gapss <- c(Hsp_gapss, Hsp_gaps)
      Hsp_align_lens <- c(Hsp_align_lens, Hsp_align_len)
      Hsp_qseqs <- c(Hsp_qseqs, Hsp_qseq)
      Hsp_hseqs <- c(Hsp_hseqs, Hsp_hseq)
      Hsp_midlines <- c(Hsp_midlines, Hsp_midline)


    }
  }

  # Create a data frame
  result_df <- data.frame(Iteration_iter_num = unlist(iter_nums),
                          Hit_num = unlist(Hit_nums),
                          Hit_id = unlist(Hit_ids),
                          Hit_def = unlist(Hit_defs),
                          Hit_accession = unlist(Hit_accessions),
                          Hit_len = unlist(Hit_lens),
                          Hsp_num = unlist(Hsp_nums),
                          Hsp_bit_score = unlist(Hsp_bit_scores),
                          Hsp_score = unlist(Hsp_scores),
                          Hsp_evalue = unlist(Hsp_evalues),
                          Hsp_query_from = unlist(Hsp_query_froms),
                          Hsp_query_to = unlist(Hsp_query_tos),
                          Hsp_hit_from = unlist(Hsp_hit_froms),
                          Hsp_hit_to = unlist(Hsp_hit_tos),
                          Hsp_query_frame = unlist(Hsp_query_frames),
                          Hsp_hit_frame = unlist(Hsp_hit_frames),
                          Hsp_identity = unlist(Hsp_identitys),
                          Hsp_positive = unlist(Hsp_positives),
                          Hsp_gaps = unlist(Hsp_gapss),
                          Hsp_align_len = unlist(Hsp_align_lens),
                          Hsp_qseq = unlist(Hsp_qseqs),
                          Hsp_hseq = unlist(Hsp_hseqs),
                          Hsp_midline = unlist(Hsp_midlines),
                          gi_id = strsplit(unlist(Hit_ids), "|", fixed = TRUE) %>% lapply(function(x) x[2]) %>% unlist())
  return(result_df)
}



check_species_match <- function(query, taxa_table) {
  for (i in seq_len(nrow(taxa_table))) {
    if (query$species == taxa_table$species[i]) {
      return(taxa_table[i,1:24])
      break
    }
  }
  return(data.frame())
}


check_species_match2 <- function(query, taxa_table, id) {
  datalist <- list()
  n = 1
  for (i in seq_len(max(taxa_table$Iteration_iter_num))) {
    current_df <- taxa_table %>% filter(Iteration_iter_num == i)
    for (j in seq_len(nrow(current_df))) {
      if (query$species == current_df$species[j]) {
        querylist <- as.list(current_df[j,1:24])
        querylist <- c(query_name = id, querylist)
        datalist[[n]] <- querylist
        n <- n+1
        break
      }
    }
  }
  return(as.data.frame(do.call(rbind, datalist)))
}


gi_to_ti_df <- function(df) {
  tax_id <- sapply(taxize::genbank2uid(unique(df$gi_id), db = "ncbi"), "[[",1)
  gi_dict <- cbind(as.list(unique(df$gi_id)))
  gi_dict <- as.data.frame(cbind(gi_dict, tax_id))

  df <- merge(df, gi_dict, by.x="gi_id", by.y="V1")
  return(df)
}

lowest_e_value <- function(df) {
  datalist <- list()
  for (i in seq_len(max(df$Iteration_iter_num))) {
    current_df <- df %>% filter(Iteration_iter_num == i)
    min_df <- current_df[current_df$Hsp_evalue == min(current_df$Hsp_evalue),]
    current_list <- as.list(min_df[1,])
    datalist[[i]] <- current_list
  }
  return(as.data.frame(do.call(rbind, datalist)))
}

print_result <- function(query, species_counts) {
  total_species <- sum(species_counts$n)
  matching_species <- species_counts[species_counts$species == query, 2] / total_species
  human_reads <- species_counts[species_counts$species == "Homo sapiens", 2] / total_species
  list((paste(query, "matching reads percent = ", as.numeric(matching_species))),
       paste(query, "human reads percent = ", as.numeric(human_reads)))
}



blast_single_result <- function(results_table, bam_file, which_result = 1, num_reads = 10, hit_list = "1") {
  final_df <- NULL
  res <- tryCatch( #If any errors, should just skip the organism
    {
      id <- results_table[which_result,1]
      if (!quiet) message("Current id: ", id)
      ti <- strsplit(id, "|", fixed = TRUE)[[1]][2]
      if (!quiet) message("Current ti: ", ti)
      fasta_seqs <- getSeqs2(id = id, bamFile = bam_file, n = num_reads)
      if (!quiet) message("Current fasta_seqs:  ", fasta_seqs)
      df <- blastXmlDf(blastSeq(x = fasta_seqs, hitListSize = as.character(hit_list), attempts = 40))
      df <- gi_to_ti_df(df)
      df
    },
    error = function(e) {
      cat("Error", conditionMessage(e))
    }
  )
  print(res)
  if (is.null(final_df)) {
    final_df <- res
  } else {
    final_df <- rbind(final_df, res)
  }

  return(final_df)
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
                    custom_format ="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \t stitle",
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


#' Blast Result Metrics
#'
#' @param blast_results_table_path
#'
#' @return Loads blast result table (csv) and generates a list of best_hit, uniqueness_score,
#' percentage_hit and contaminant score
#'
#' TODO: ADD ERROR HANDLING

blast_result_metrics <- function(blast_results_table_path,
                                 uniqueness_score_by = "species",
                                 percentage_hit_by = "species",
                                 contaminant_score_by = "genus"){
  tryCatch(
    {
      # Load in blast results table
      blast_results_table <- read.csv(blast_results_table_path, header = TRUE)

      # Remove any empty tables
      if (nrow(blast_results_table) == 0) {
        return(list(best_hit = 0, uniqueness_score = 0, percentage_hit = 0, contaminant_score = 0))
      }

      # Get the best blast result for every query (read)
      lowest_eval_per_query <- blast_results_table %>% group_by(qseqid) %>%
        slice_min(evalue, with_ties = FALSE)
      summary_table <- blast_results_table %>% group_by(name) %>% summarise(num_reads = n())
      summary_table$genus <- strsplit(summary_table$name, split = " ")[[1]][1]
      summary_table$species <- strsplit(summary_table$name, split = " ")[[1]][2]
      summary_table_genus <- summary_table %>% group_by(genus) %>% summarise(num_reads = sum(num_reads))
      summary_table_species <- summary_table %>% group_by(species) %>% summarise(num_reads = sum(num_reads))

      best_hit <- summary_table %>% slice_max(num_reads, with_ties = FALSE)
      assign(paste0("best_hit_", percentage_hit_by),
             eval(parse(text = paste0("summary_table_", percentage_hit_by))) %>%
               slice_max(num_reads, with_ties = FALSE))

      # Calculate Metrics
      uniqueness_score <- 1/nrow(eval(parse(text = paste0("summary_table_", uniqueness_score_by))))

      percentage_hit <- eval(parse(text = paste0("best_hit_", percentage_hit_by)))$num_reads /
        sum(eval(parse(text = paste0("summary_table_", percentage_hit_by)))$num_reads)

      contaminant_score <- eval(parse(text = paste0("summary_table_", contaminant_score_by))) %>%
        arrange(desc(num_reads)) %>% slice(2) %>%
        select(2) / sum(eval(parse(text = paste0("summary_table_", contaminant_score_by)))
                        %>% select(num_reads))
      return(list(best_hit = best_hit$name, uniqueness_score = uniqueness_score, percentage_hit = percentage_hit, contaminant_score = unlist(contaminant_score)))
    },
    error = function(e)
    {
      cat("Error", conditionMessage(e), "/n")
      return(list(best_hit = 0, uniqueness_score = 0, percentage_hit = 0, contaminant_score = 0))
    }
  )
}








