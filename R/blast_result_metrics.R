#' Blast Result Metrics
#'
#' @param blast_results_table_path path for blast results csv file
#'
#' @return Loads blast result table (csv) and generates a list of best_hit, uniqueness_score,
#' percentage_hit and contaminant score
#'

blast_result_metrics <- function(blast_results_table_path){
  tryCatch(
    {
      # Load in blast results table
      blast_results_table <- read.csv(blast_results_table_path, header = TRUE)

      # Remove any empty tables
      if (nrow(blast_results_table) < 2) {
        return(data.frame(best_hit = 0,
                          uniqueness_score = 0,
                          species_percentage_hit = 0,
                          genus_percentage_hit = 0,
                          species_contaminant_score = 0,
                          genus_contaminant_score = 0))
      }

      # Adding Species and Genus columns
      blast_results_table <- blast_results_table %>%
        mutate(MetaScope_genus = word(MetaScope_Genome, 1, 1, sep = " ")) %>%
        mutate(MetaScope_species = word(MetaScope_Genome, 1, 2, sep = " ")) %>%
        mutate(query_genus = word(name, 1, 1, sep = " ")) %>%
        mutate(query_species = word(name, 1, 2, sep = " "))

      # Removing duplicate query num and query species
      blast_results_table <- blast_results_table %>%
        distinct(qseqid, query_species, .keep_all = TRUE)

      # Calculating Metrics
      best_hit <- blast_results_table %>%
        group_by(query_species) %>%
        summarise(num_reads = n()) %>%
        slice_max(num_reads, with_ties = FALSE)

      uniqueness_score <- blast_results_table %>%
        group_by(query_species) %>%
        summarise(num_reads = n()) %>%
        nrow()

      species_percentage_hit <- blast_results_table %>%
        filter(MetaScope_species == query_species) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      genus_percentage_hit <- blast_results_table %>%
        filter(MetaScope_genus == query_genus) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      species_contaminant_score <- blast_results_table %>%
        filter(MetaScope_species != query_species) %>%
        distinct(qseqid, .keep_all = TRUE) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      genus_contaminant_score <- blast_results_table %>%
        filter(MetaScope_genus != query_genus) %>%
        distinct(qseqid, .keep_all = TRUE) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      return(data.frame(best_hit = best_hit$query_species,
                        uniqueness_score = uniqueness_score,
                        species_percentage_hit = species_percentage_hit,
                        genus_percentage_hit = genus_percentage_hit,
                        species_contaminant_score = species_contaminant_score,
                        genus_contaminant_score = genus_contaminant_score))
    },
    error = function(e)
    {
      cat("Error", conditionMessage(e), "/n")
      return(data.frame(best_hit = 0,
                        uniqueness_score = 0,
                        species_percentage_hit = 0,
                        genus_percentage_hit = 0,
                        species_contaminant_score = 0,
                        genus_contaminant_score = 0))
    }
  )
}
