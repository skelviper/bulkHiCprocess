  
#! /home/skelviper/anaconda3/envs/R/bin/R
# change above to your own R path
# scripts for clean compartment-type of TAD in isolation calling method.
#@author zliu
#@2021.2.1

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | args[1]=="-h") {
  stop("usage: Rscript cleanTADs.R compartment.bed standardTAD TADclean.bed", call.=FALSE)
}

library(tidyverse)

com <- read.delim(args[1],header = F,col.names = c("chrom","start1","end1","type"))
rawTADs <- read_table2(args[2])

processRawTAD <- function(file,boundary_strength_limit = 0.1){
    temp<-rawTADs %>% filter(boundary_strength_200000 > 0.1 & is_bad_bin == FALSE) %>% select(1,2)
    temp <- file %>% filter(boundary_strength_200000 > boundary_strength_limit & is_bad_bin == FALSE) %>% left_join(temp,by="chrom") %>% filter(start.y-start.x<1500000 & start.x < start.y)
    temp <- temp %>% select(-end) %>% select(chrom,start.x,start.y,everything())
    return(temp)
}

tad <- processRawTAD(rawTADs) %>% select(1,2,3)
names(tad) <- c("chrom","start","end")

tadclean <- setdiff(tad,tad %>% full_join(com) %>% arrange(desc(type)) %>% filter((abs(start-start1)<=100000&abs(end-end1)<=100000)) %>% select(chrom,start,end))

write_tsv(tadclean,args[3])