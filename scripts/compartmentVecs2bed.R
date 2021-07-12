  
#! /home/skelviper/anaconda3/envs/R/bin/R
# change above to your own R path
# scripts for convert cooltools compartment vecs to bed file.
#@author zliu
#@2021.2.1

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | args[1]=="-h") {
  stop("usage: Rscript compartmentVecs2bed.R compartment.vecs compartment.bed", call.=FALSE)
}

library(tidyverse)

compartmentVecs <- read_table2(args[1])
com <-compartmentVecs %>% select(chrom,start,end,E1) %>% mutate(compartment = ifelse(E1>0,"A","B"),strand_as_com = ifelse(E1>0,"+","-") )%>% na.omit()

write_tsv(com,args[2],col_names = F)