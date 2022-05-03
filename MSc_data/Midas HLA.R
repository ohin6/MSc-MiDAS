knitr::opts_chunk$set(collapse = TRUE)

# load packages
library("tidyverse")
library("kableExtra")
library("knitr")
library("midasHLA")

# set working directory
setwd("/Users/owen/Library/Mobile Documents/com~apple~CloudDocs/midasHLA/MSc_data/")

# import files
map = tibble(read.table("HLA_files/hla-subset.map"))

genotypes = tibble(read.table("HLA_files/hla-subset.ped")) %>%
  select(-c(V2:V6))

#####change colnames in genotype file #####

#Get row names
colnames = map$V2
col_append = "ID"

# loop function to get diploid version of each genotype
for (i in 1:length(colnames))
{
  x = paste(colnames[i],"a")
  y = paste(colnames[i],"b")
  col_append = append(col_append, x)
  col_append = append(col_append, y)
}
col_append = gsub(" ", "", col_append, fixed = TRUE) # remove white spaces

colnames(genotypes) = col_append

#### convert into MiDAs format####
