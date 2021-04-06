library(glue)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(splitstackshape)
library(readr)
library(data.table)
library(rebus)

name_key <- read.delim("~/projects/b1059/projects/Ryan/csq/flat_file/wormbase_name_key.txt") %>%
  select(WormBase.Gene.ID, Public.Name)

data<- data.table::fread("~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile.tsv")


add_gene <- data %>%
  left_join(name_key, by = c( "GENE" = "WormBase.Gene.ID")) %>%
  rename("WORMBASE_ID" = "GENE") %>% rename("GENE" = "Public.Name")


data.table::fwrite(add_gene, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile-gene.tsv")
