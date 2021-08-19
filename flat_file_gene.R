library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(data.table)


name_key <- data.table::fread("~/projects/b1059/projects/Ryan/csq/flat_file/wormbase_name_key.txt") %>%
  dplyr::select(WormBase.Gene.ID, Public.Name)

data<- data.table::fread("~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile.tsv")


add_gene <- data %>%
  dplyr::left_join(name_key, by = c( "GENE" = "WormBase.Gene.ID")) %>%
  dplyr::rename("WORMBASE_ID" = "GENE") %>% rename("GENE" = "Public.Name")


data.table::fwrite(add_gene, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile-gene.tsv")
