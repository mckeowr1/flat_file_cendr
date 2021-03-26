library(dplyr)
library(tidyr)
library(data.table)
library(readr)

#Reformat table
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

sample_BCSQ <- data.table::fread("~/projects/b1059/projects/Sophie/csq/re-an_strain_BCSQ.tsv")

sample_parsed <- parse_strain(sample_BCSQ)

#Add scores
score_unparsed_annotation <- read.delim("~/projects/b1059/projects/Sophie/csq/score_unparsed_annotation.tsv")

strain_variant_score <- dplyr::left_join(sample_parsed, score_unparsed_annotation)

strain_variant_score <- strain_variant_score %>%
  tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|") %>%
  na_if("NA")


#Remove non-single AA substitutions. Make columns integers
clean_flat_file <- function(df){
  clean_flat_file <- df %>%
    dplyr::mutate(BLOSUM = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" & 
                                    CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, BLOSUM)) %>%
    dplyr::mutate(Grantham = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" & 
                                      CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, Grantham)) %>%
    dplyr::mutate(Grantham = sapply(clean_flat_file$Grantham, as.integer)) %>% 
    dplyr::mutate(BLOSUM = sapply(clean_flat_file$BLOSUM, as.integer)) %>%
    dplyr::mutate(Percent_Protein = sapply(clean_flat_file$Percent_Protein, as.integer))
}

cleaned_flat_file <- clean_flat_file(strain_variant_score)
  

readr::write_tsv(cleaned_flat_file, "~/projects/b1059/projects/Sophie/csq/strain_variant_table.tsv")

