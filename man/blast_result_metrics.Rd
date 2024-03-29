% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast_result_metrics.R
\name{blast_result_metrics}
\alias{blast_result_metrics}
\title{Blast Result Metrics}
\usage{
blast_result_metrics(
  blast_results_table_path,
  uniqueness_score_by = "species",
  percentage_hit_by = "species",
  contaminant_score_by = "genus"
)
}
\arguments{
\item{blast_results_table_path}{path for blast results csv file}

\item{uniqueness_score_by}{score between 0-1 for diversity of blast hits
\itemize{
\item "species" (the default):  scores uniqueness by species
\item "genus": scores uniqueness by genus
}}

\item{contaminant_score_by}{percentage of blast results that matches the second best hit
\itemize{
\item "species" (the default):  scores percentage hit by species
\item "genus": scores percentage hit by genus
}}

\item{percetange_hit_by}{percentage of blast results that matches best hit
\itemize{
\item "species" (the default):  scores percentage hit by species
\item "genus": scores percentage hit by genus
}}
}
\value{
Loads blast result table (csv) and generates a list of best_hit, uniqueness_score,
percentage_hit and contaminant score
}
\description{
Blast Result Metrics
}
