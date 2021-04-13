library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(rebus)




order_annotation <- function(df){
  df[with(df, order(CHROM, POS)), ]
}


impact_scoring <- function(df) {  #Score 3 or Less
  multi_con<- df %>% #Subset multi-consequence variants
  dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&"))) %>%
    dplyr::filter(.$multi_con == TRUE) %>% select(!multi_con) %>% as_tibble()

  refromated_multi <- multi_con %>%
    tidyr::separate(CONSEQUENCE, sep = "&", into = c("CON1", "CON2", "CON3"), remove= FALSE) %>% #Separate multi consequences
    dplyr::mutate(CON1 = sapply(.$CON1, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    dplyr::mutate(CON2 = sapply(.$CON2, impact_numeric)) %>% #Recode impact to numeric 4 - highest 1- lowest
    dplyr::mutate(CON3 = sapply(.$CON3, impact_numeric)) %>%#Recode impact to numeric 4 - highest 1- lowest
    dplyr::mutate(VARIANT_IMPACT = pmax(.$CON1, .$CON2, .$CON3, na.rm = TRUE)) %>% #Select the highest impact to be variant impact
    dplyr::mutate(VARIANT_IMPACT = sapply(VARIANT_IMPACT, impact_numeric_tocharacter)) %>% #Convert from numeric to character
    dplyr::mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character)) %>%
    select(!CON1, !CON2)%>% #Remove columns created in previous steps
    as.data.table()

  single_con <- df %>%
    dplyr::mutate(multi_con = stringr::str_detect(.$CONSEQUENCE, pattern = fixed("&"))) %>% #Mark multi-consequence variants
    dplyr::filter(.$multi_con == FALSE) %>% #Subset single consequence
    dplyr::mutate(VARIANT_IMPACT = sapply(.$CONSEQUENCE, impact)) %>% #Add impact score
    dplyr::select(!multi_con) %>% #Get rid of multi-consequence flag
    dplyr::mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character)) %>%
    as.data.table()

  clean <- data.table::rbindlist(list(single_con,refromated_multi), fill=TRUE) %>% #Bind multi_con and single con back together
    order_annotation() %>% #Re-order
    dplyr::mutate(VARIANT_IMPACT = sapply(.$VARIANT_IMPACT, as.character))
  return(clean)
}
#Functions called above
impact <- function(x) {  #Converts consequence impact valye
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
  dplyr::recode(x,
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

}
impact_numeric_tocharacter <- function(x){
  dplyr::recode(x,
         "4"="HIGH",
         '3'="HIGH",
         "2"="LOW",
         '1'="LOW")
}


data <-data.table::fread("~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile-gene.tsv")



impacted <- impact_scoring(data) %>% #Perfrom Impact Scoring
  dplyr::mutate(linker = stringr::str_detect(.$CONSEQUENCE, pattern = "@" )) #Linker check for Linker impact score


impacted$VARIANT_IMPACT = if_else(impacted$linker == TRUE, "Linker", impacted$VARIANT_IMPACT ) #Add linker variant impact



table_ready <- impacted %>%
  dplyr::select(-c(linker, CON1, CON2, CON3))

data.table::fwrite(table_ready, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile-gene-impact.tsv")
