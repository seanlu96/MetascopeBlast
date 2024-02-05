#' Blast Result Metrics
#'
#' @param blast_results_table_path path for blast results csv file
#' @param uniqueness_score_by score between 0-1 for diversity of blast hits
#'  * "species" (the default):  scores uniqueness by species
#'  * "genus": scores uniqueness by genus
#'@param percetange_hit_by percentage of blast results that matches best hit
#'  * "species" (the default):  scores percentage hit by species
#'  * "genus": scores percentage hit by genus
#'@param contaminant_score_by percentage of blast results that matches the second best hit
#'  * "species" (the default):  scores percentage hit by species
#'  * "genus": scores percentage hit by genus
#'
#' @return Loads blast result table (csv) and generates a list of best_hit, uniqueness_score,
#' percentage_hit and contaminant score
#'

blast_result_metrics <- function(blast_results_table_path,
                                 uniqueness_score_by = "species",
                                 percentage_hit_by = "species",
                                 contaminant_score_by = "genus"){
  tryCatch(
    {
      # Load in blast results table
      blast_results_table <- read.csv(blast_results_table_path, header = TRUE)

      # Remove any empty tables
      if (nrow(blast_results_table) < 2) {
        return(data.frame(best_hit = 0, uniqueness_score = 0, percentage_hit = 0, contaminant_score = 0))
      }

      # Get the best blast result for every query (read)
      lowest_eval_per_query <- blast_results_table %>% group_by(qseqid) %>%
        slice_min(evalue, with_ties = FALSE)
      summary_table <- blast_results_table %>% group_by(name) %>% summarise(num_reads = n())
      summary_table$genus <- strsplit(summary_table$name, split = " ")[[1]][1]
      summary_table$species <- strsplit(summary_table$name, split = " ")[[1]][1:2]
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

      if (eval(parse(text = paste0("summary_table_", contaminant_score_by))) %>% nrow() == 1) {
        contaminant_score <- 0
      } else{
        contaminant_score <- eval(parse(text = paste0("summary_table_", contaminant_score_by))) %>%
          arrange(desc(num_reads)) %>% slice(2) %>%
          select(2) / sum(eval(parse(text = paste0("summary_table_", contaminant_score_by)))
                          %>% select(num_reads))
      }


      return(data.frame(best_hit = best_hit$MetaScope_Genome,
                        uniqueness_score = uniqueness_score,
                        percentage_hit = percentage_hit,
                        contaminant_score = contaminant_score))
    },
    error = function(e)
    {
      cat("Error", conditionMessage(e), "/n")
      return(data.frame(best_hit = 0, uniqueness_score = 0, percentage_hit = 0, contaminant_score = 0))
    }
  )
}
