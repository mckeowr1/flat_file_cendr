library(tidyverse)
library(dplyr)
library(tidyr)
library(magrittr)
library(valr)
library(fuzzyjoin)
library(Bioconductor)

#Load Divergent Regions
divergent <- data.table::fread("lee2020.divergent_regions_strain.bed", 
                               col.names = c("chrom", "start", "end", "strain"))

#Condense Overlapping Regions (If portion of regions shared by >1 strian 
condensed <- divergent %>%
  valr::bed_merge() %>% 
  dplyr::mutate(DIVERGENT = "D") # Add marker to Divergent Region 

#Read In Flat File 
data <- data.table::fread("WI.20210121.strain.annotation-GB-gene-impact.tsv")


#Join the data - if a position is within a divergent region Divergent tag is added

join <- genome_join(data, condensed,  by = c(
  "CHROM" = "chrom",
  "POS" = "start",
  "POS" = "end"), 
  mode = "left") %>%  
  select(-chrom, -start, -end)




data.table::fwrite(join, "WI.20210121.strain.annotation-GB-gene-impact-divergent.tsv" )
