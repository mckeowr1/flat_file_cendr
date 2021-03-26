library(glue)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(splitstackshape)
library(readr)
library(data.table)
library(rebus)




order_annotation <- function(df){
  df[with(df, order(CHROM, POS)), ]
}


# parse_VCF <- function(df){
#   parsed <- df %>%
#     dplyr::rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4", "SAMPLE" = "V5", "ANNOTATION" = "V6") %>%
#     tidyr::separate_rows(ANNOTATION, sep=",")%>%
#     tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|")
# }

impact_scoring <- function(df) {  #Works for up to two multi consequences - can handle 3 but the highest may not always be resulted
  multi_con<- df %>% dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&"))) %>%
    filter(.$multi_con == TRUE) %>% select(!multi_con) %>% as_tibble()  #Filter for only multiconsequnce variants
  
  refromated_multi <- multi_con %>%
    separate(CONSEQUENCE, sep = "&", into = c("CON1", "CON2", "CON3"), remove= FALSE) %>% #Separate consequnces (Works for 2 ATM)
    mutate(CON1 = sapply(.$CON1, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    mutate(CON2 = sapply(.$CON2, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    mutate(CON3 = sapply(.$CON3, impact_numeric)) %>%
    mutate(VARIANT_IMPACT = pmax(.$CON1, .$CON2, .$CON3, na.rm = TRUE)) %>% #Select the highest impact to be variant impact
    mutate(VARIANT_IMPACT = sapply(VARIANT_IMPACT, impact_numeric_tocharacter)) %>% #Convert from numeric to character
    mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character)) %>% select(!CON1, !CON2)%>% as.data.table()
  
  single_con <- df %>% dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&"))) %>%
    filter(.$multi_con == FALSE) %>% dplyr::mutate(VARIANT_IMPACT = sapply(.$CONSEQUENCE, impact)) %>%
    select(!multi_con) %>% mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character)) %>% as.data.table()
  
  clean <- data.table::rbindlist(list(single_con,refromated_multi), fill=TRUE) %>%
    order_annotation() %>% mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character))
  return(clean)
} #Uses Below functions to handle impact conversion & Multi-consequence variants
impact <- function(x) {  #Edited to group Moderat and High Into High
  dplyr::recode(x,
                
                "missense" ="HIGH",
                "synonymous"="LOW",
                "stop_lost"= "HIGH",
                "stop_gained"="HIGH",
                "inframe_deletion"="HIGH",
                "inframe_insertion"="HIGH",
                "frameshift"="HIGH",
                "splice_acceptor"="HIGH",
                "splice_donor"="HIGH",
                "start_lost"="HIGH",
                "splice_region"="HIGH",
                "stop_retained"="LOW",
                "5_prime_utr"="LOW",
                "3_prime_utr"="LOW",
                "non_coding"="LOW",
                "intron"="LOW",
                "intergenic"="LOW",
                "inframe_altering"="HIGH",
                "coding_sequence"="LOW",
                "feature_elongation"="LOW",
                "start_retained"="LOW" ,
                "*missense" ="HIGH",
                "*synonymous"="LOW",
                "*stop_lost"= "HIGH",
                "*stop_gained"="HIGH",
                "*inframe_deletion"="HIGH",
                "*inframe_insertion"="HIGH",
                "*frameshift"="HIGH",
                "*splice_acceptor"="HIGH",
                "*splice_donor"="HIGH",
                "*start_lost"="HIGH",
                "*splice_region"="HIGH",
                "*stop_retained"="LOW",
                "*5_prime_utr"="LOW",
                "*3_prime_utr"="LOW",
                "*non_coding"="LOW",
                "*intron"="LOW",
                "*intergenic"="LOW",
                "*inframe_altering"="HIGH",
                "*coding_sequence"="LOW",
                "*feature_elongation"="LOW",
                "*start_retained"="LOW")
  
} #Converts Consequence to Impact
impact_numeric <- function(x) {
  recode(x,
         "missense" =3,
         "synonymous"=1,
         "stop_lost"= 4,
         "stop_gained"=4,
         "inframe_deletion"=3,
         "inframe_insertion"=3,
         "frameshift"=4,
         "splice_acceptor"=4,
         "splice_donor"=4,
         "start_lost"=4,
         "splice_region"=3,
         "stop_retained"=1,
         "5_prime_utr"=2,
         "3_prime_utr"=2,
         "non_coding"=2,
         "intron"=2,
         "intergenic"=2,
         "inframe_altering"=3,
         "coding_sequence"=2,
         "feature_elongation"=2,
         "start_retained"=1 ,
         "*missense" =3,
         "*synonymous"=1,
         "*stop_lost"= 4,
         "*stop_gained"=4,
         "*inframe_deletion"=3,
         "*inframe_insertion"=3,
         "*frameshift"=4,
         "*splice_acceptor"=4,
         "*splice_donor"=4,
         "*start_lost"=4,
         "*splice_region"=3,
         "*stop_retained"=1,
         "*5_prime_utr"=2,
         "*3_prime_utr"=2,
         "*non_coding"=2,
         "*intron"=2,
         "*intergenic"=2,
         "*inframe_altering"=3,
         "*coding_sequence"=2,
         "*feature_elongation"=2,
         "*start_retained"=1)
  
}  #Handles multi-consequence with numneric to character
impact_numeric_tocharacter <- function(x){
  recode(x,
         "4"="HIGH",
         '3'="HIGH",
         "2"="LOW",
         '1'="LOW")
}


data <-data.table::fread("WI.20210121.strain.annotation-GB-gene.tsv")



impacted <- impact_scoring(data) %>% #Perfrom Impact Scoring
  mutate(linker = str_detect(.$CONSEQUENCE, pattern = "@" )) #Linker check for Linker impact score

print("General Imapct")

impacted$VARIANT_IMPACT = if_else(impacted$linker == TRUE, "Linker", impacted$VARIANT_IMPACT )

print("Linkers Impact")

table_ready <- impacted %>%  
  select(-c(linker, CON1, CON2, CON3))

data.table::fwrite(table_ready, "WI.20210121.strain.annotation-GB-gene-impact.tsv")




