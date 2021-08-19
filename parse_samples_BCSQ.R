library(tidyverse)
library(data.table)

#Requires at least 32 GB mem in slurm script runs for >2hrs

parse_strain <- function(df){
  parsed_strain <- df %>%
    dplyr::rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4") %>%
    dplyr::mutate(V5=gsub('.{1}$', '', df$V5)) %>%
    tidyr::separate_rows(V5, sep = "=")%>%
    tidyr::separate(V5, into = c("Strain", "ANNOTATION"), sep = ":")%>%
    tidyr::separate_rows(ANNOTATION, sep = ",") %>%
    dplyr::group_by(CHROM, POS, REF, ALT, ANNOTATION) %>%
    dplyr::summarise(Strains = paste(Strain, collapse = ","))
}

sample_BCSQ <- data.table::fread("~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401_strain_bcsq.tsv")

sample_parsed <- parse_strain(sample_BCSQ)

readr::write_tsv(sample_parsed, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile.tsv")