# Coverage Plot
# Brie Odom
# 8/23/21

# Function to check the locations of various species

# Setup -----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Rsamtools)
  library(tidyverse)
  library(ggpubr)
})

# Helper function -------------------------------------------------------------
locations <- function(taxid, input_bam) {
  df <- as.data.frame(Rsamtools::seqinfo(input_bam))
  df$seqnames <- row.names(df)
  allGenomes <- grep(paste0("ti|", taxid), df$seqnames, value = TRUE, fixed = TRUE)
  Genome <- sample(allGenomes, 1)
  print(Genome)
  species_name <- strsplit(Genome, split = "|", fixed = TRUE)[[1]][4]
  p2 <- ScanBamParam(what = c("rname", "pos", "qname", "qwidth", "seq"),
                     which = GRanges(Genome, IRanges(1, 1e+07)))

  # "Upload" bam file into environment
  res0 <- scanBam(input_bam, param = p2)[[1]]

  p <- ggplot(mapping = aes(x=res0$pos)) +
    geom_histogram(binwidth = 1000) +
    labs(x = "Leftmost position in genome",
        title = species_name,
        y = "read count")

  aa <- quantile(res0$pos)
  tab <- aa %>% as.data.frame() %>% t() %>% ggtexttable()
  re <- list(plot = p, table = tab)
  return(re)
}

# Main function ---------------------------------------------------------------
get_locations <- function(input_bam, species) {
  # Obtain plots, quantiles of locations
  re <- lapply(species,
               function(x) locations(x, input_bam))
  names(re) <- species
  return(re)
}

# Example ---------------------------------------------------------------------
input_bam = Rsamtools::BamFile("E:/bls11_S46_sorted.bam", index = "E:/bls11_S46_sorted.bam")
test_species = c(652616, 419947, 1031709, 83331, 1097669, 443149)
aa=get_locations(input_bam, # Location of bam file
              test_species)

ggpubr::ggarrange(aa[[1]]$plot,
                  aa[[2]]$plot,
                  aa[[3]]$plot,
                  aa[[4]]$plot,
                  aa[[5]]$plot,
                  aa[[6]]$plot,
                  nrow = 3, ncol = 2)
