###### Before starting script make sure you change working directory #######


knitr::opts_chunk$set(collapse = TRUE)

## Downlaod or load relevant pacakages

BiocManager::install("midasHLA")


library("tidyverse")
library("kableExtra")
library("knitr")
library("midasHLA")
library('stringr')

## set working directory
#setwd("/Users/owen/Library/Mobile Documents/com~apple~CloudDocs/midasHLA/MSc_data/")


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

#transpose table


TransGenotypes = as.data.frame(t(genotypes))

colnames(TransGenotypes) = genotypes$ID
TransGenotypes = TransGenotypes[-1,]

# filter rows containing P

patientIds = genotypes$ID

# create empty tibble which contains different HLA genes
MiDASGeno = tibble("HLA_A1", "HLA_A2","HLA_A3","HLA_A4", "HLA_C1", "HLA_C2","HLA_C3","HLA_C4", "HLA_B1", "HLA_B2","HLA_B3","HLA_B4", "HLA_DRB1", "HLA_DRB12","HLA_DRB13","HLA_DRB14","HLA_DQA11", "HLA_DQA12","HLA_DQA13","HLA_DQA14", "HLA_DQB11", "HLA_DQB12","HLA_DQB13","HLA_DQB14", "HLA_DPA11", "HLA_DPA12","HLA_DPA13","HLA_DPA14","HLA_DPB11", "HLA_DPB12","HLA_DPB13","HLA_DPB14")
MiDASGeno = MiDASGeno[-1,]
dfRowNames = vector()

# For loop through each column in Transgenotype dataframe which contains 'p' for present allele


for (i in seq_along(patientIds)){
  df = TransGenotypes[i] %>%
    filter(str_detect(TransGenotypes[[i]], "P"))
  dfRowNames = row.names(df)
  MiDASGeno = rbind(MiDASGeno,dfRowNames)
}

MiDASGeno2 = MiDASGeno
rownames(MiDASGeno2) = patientIds

#### tidying data ###

#remove last letter from every element
for (i in 1:ncol(MiDASGeno2)){
  MiDASGeno2[[i]] = gsub("a$","", MiDASGeno2[[i]])
  MiDASGeno2[[i]] = gsub("b$","", MiDASGeno2[[i]])
}

#transpose
MiDASGeno3 = t(MiDASGeno2)
MiDASGeno3 = as.data.frame(MiDASGeno3)

###### need to remove 2 digit HLA######

MiDASGeno4=tibble()

for (i in 1:patientIds){
  HLA = MiDASGeno3[i] %>%
    filter(str_detect(MiDASGeno3[[i]], "_\\d{4}$"))
  rowbind = HLA[[1]]
  MiDASGeno4 = rbind(MiDASGeno4,rowbind)
}
colnames(MiDASGeno4) = c("A_1", "A_2", "C_1", "C_2", "B_1", "B_2", "DRB1_1", "DRB1_2","DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", "DPA1_1", "DPA1_2","DPB1_1", "DPB1_2")


##### Because they're seems to be errors genotyping file where numerous alleles are detected, a bit of fidging needs to be done






MiDASGeno4 = MiDASGeno4 %>%
  select(A_1,	A_2,	B_1,	B_2,	C_1,	C_2,	DPA1_1,	DPA1_2,	DPB1_1,	DPB1_2,	DQA1_1,	DQA1_2,	DQB1_1,	DQB1_2,	DRB1_1,	DRB1_2)
  
df = MiDASGeno4 %>%
  mutate(across(everything(), str_replace_all, 'HLA_', '')) %>%
  mutate(across(everything(), str_replace_all, '_', '*'))
  

### Create function to add character
fun_insert = function(x, pos, insert){
  gsub(paste0("^(.{", pos,"})(.*)$"),
       paste0("\\1",insert,"\\2"),
       x)
}
df=tibble(df)

df = mutate(df, across(everything(), ~ fun_insert(.x,str_length(.x)-2,':')))



write_tsv(MiDASGeno4, "MiDASGeno.txt")


map %>%
filter(str_detect(map$V2, 'HLA_DRB'))
