library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(hablar)

#Read In Flat File

sample_parsed <- data.table::fread("WI.20210121.strain.annotation-integer.tsv")        # parse_strain(sample_BCSQ)
  
#Read unparsed scored annotations 
# score_unparsed_annotation <- read.delim("score_unparsed_annotation.tsv")
# 
# strain_variant_score <- dplyr::left_join(sample_parsed, score_unparsed_annotation)
# 
# strain_variant_score <- strain_variant_score %>%
#   tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|") %>%
#   na_if("NA")

clean_flat_file <- function(df){
  clean_flat_file <- df %>%
    dplyr::mutate(BLOSUM = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" & 
                                    CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, BLOSUM)) %>%
    dplyr::mutate(Grantham = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" & 
                                      CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, Grantham)) %>%
    hablar::convert(num(BLOSUM, Grantham, Percent_Protein))
    # dplyr::mutate(Grantham = sapply(clean_flat_file$Grantham, as.integer)) %>% 
    # dplyr::mutate(BLOSUM = sapply(clean_flat_file$BLOSUM, as.integer)) %>%
    # dplyr::mutate(Percent_Protein = sapply(clean_flat_file$Percent_Protein, as.integer))
}


cleaned_flat_file <- clean_flat_file(sample_parsed)



readr::write_tsv(cleaned_flat_file, "WI.20210121.strain.annotation-GB.tsv")
