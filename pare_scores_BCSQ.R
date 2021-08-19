library(dplyr)
library(tidyr)
library(data.table)
library(readr)

scores <- data.table::fread("/projects/b1059/projects/Sophie/DataExplorers/csq/MakingFlatFileCode/Score-Annotation.tsv")

parse_score <- function(df){
  parsed <- df %>%
    dplyr::rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4", "ANNOTATION" = "V5", "BLOSUM" = "V6", "Grantham" = "V7", "Percent_Protein" = "V8") %>%
    tidyr::separate_rows(ANNOTATION, BLOSUM, Grantham, Percent_Protein, sep=",")

}

scored_with_annotation <- parse_score(scores)

readr::write_tsv(scored_with_annotation, "/projects/b1059/projects/Sophie/DataExplorers/csq/WI.20210121.Scores.Annotation.tsv")
