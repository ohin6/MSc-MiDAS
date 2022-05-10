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


MiDASGeno4[53,] = c("HLA_A_0101", "HLA_A_0102", "HLA_C_0501", "HLA_C_1202", "HLA_B_4402", "HLA_B_5201", "HLA_DRB1_1301","HLA_DRB1_1502","HLA_DQA1_0103","HLA_DQA1_0103", "HLA_DQB1_0601","HLA_DQB1_0603","HLA_DPA1_0103","HLA_DPA1_0201","HLA_DPB1_0401","HLA_DPB1_1701")
MiDASGeno4[61,] = c('HLA_A_0201' , 'HLA_A_2402', 'HLA_C_0304 ' ,  'HLA_C_0305' ,  'HLA_B_1401' , 'HLA_B_4001', 'HLA_DRB1_0701' ,  'HLA_DRB1_0701' ,  'HLA_DQA1_0201' ,  'HLA_DQA1_0201' ,  'HLA_DQB1_0202' ,  'HLA_DQB1_0202' ,  'HLA_DPA1_0103' ,  'HLA_DPA1_0103' ,  'HLA_DPB1_0201' ,  'HLA_DPB1_0201')
MiDASGeno4[217,] = c('HLA_A_0201', 'HLA_A_2601', 'HLA_C_0702', 'HLA_C_0702', 'HLA_B_0702', 'HLA_B_4102', 'HLA_DRB1_0101', 'HLA_DRB1_1303', 'HLA_DQA1_0101', 'HLA_DQA1_0501', 'HLA_DQB1_0201', 'HLA_DQB1_0501', 'HLA_DPA1_0103', 'HLA_DPA1_0201', 'HLA_DPB1_0101', 'HLA_DPB1_0201')
MiDASGeno4[399,] = c('HLA_A_0101', 'HLA_A_1101', 'HLA_C_0501', 'HLA_C_1202', 'HLA_B_4402', 'HLA_B_5201', 'HLA_DRB1_1301', 'HLA_DRB1_1502', 'HLA_DQA1_0103', 'HLA_DQA1_0103', 'HLA_DQB1_0601', 'HLA_DQB1_0603', 'HLA_DPA1_0103', 'HLA_DPA1_0201', 'HLA_DPB1_0401', 'HLA_DPB1_1701')
MiDASGeno4[696,] = c('HLA_A_0201', 'HLA_A_1101', 'HLA_C_0701', 'HLA_C_0702', 'HLA_B_1801', 'HLA_B_1801', 'HLA_DRB1_0401', 'HLA_DRB1_0404', 'HLA_DQA1_0301', 'HLA_DQA1_0301', 'HLA_DQB1_0302', 'HLA_DQB1_0302', 'HLA_DPA1_0103', 'HLA_DPA1_0201', 'HLA_DPB1_1101', 'HLA_DPB1_2301')
MiDASGeno4[714,] = c('HLA_A_0101', 'HLA_A_1101', 'HLA_C_0202', 'HLA_C_0701', 'HLA_B_0702', 'HLA_B_0801', 'HLA_DRB1_0301', 'HLA_DRB1_0701', 'HLA_DQA1_0201', 'HLA_DQA1_0501', 'HLA_DQB1_0201', 'HLA_DQB1_0303', 'HLA_DPA1_0103', 'HLA_DPA1_0201', 'HLA_DPB1_0101', 'HLA_DPB1_0401')
MiDASGeno4[849,] = c('HLA_A_0201', 'HLA_A_2902', 'HLA_C_0602', 'HLA_C_1601', 'HLA_B_4404', 'HLA_B_5701', 'HLA_DRB1_0701', 'HLA_DRB1_0701', 'HLA_DQA1_0201', 'HLA_DQA1_0201', 'HLA_DQB1_0202', 'HLA_DQB1_0303', 'HLA_DPA1_0103', 'HLA_DPA1_0103', 'HLA_DPB1_0301', 'HLA_DPB1_0401')
MiDASGeno4[970,] = c('HLA_A_0301', 'HLA_A_6802', 'HLA_C_0401', 'HLA_C_0702', 'HLA_B_0702', 'HLA_B_5301', 'HLA_DRB1_0401', 'HLA_DRB1_1302', 'HLA_DQA1_0102', 'HLA_DQA1_0301', 'HLA_DQB1_0604', 'HLA_DQB1_0604', 'HLA_DPA1_0103', 'HLA_DPA1_0103', 'HLA_DPB1_0401', 'HLA_DPB1_0401')
MiDASGeno4[1275,] = c('HLA_A_0201', 'HLA_A_2402', 'HLA_C_0202', 'HLA_C_0303', 'HLA_B_1507', 'HLA_B_2705', 'HLA_DRB1_0404', 'HLA_DRB1_0404', 'HLA_DQA1_0101', 'HLA_DQA1_0301', 'HLA_DQB1_0302', 'HLA_DQB1_0501', 'HLA_DPA1_0103', 'HLA_DPA1_0201', 'HLA_DPB1_0601', 'HLA_DPB1_0901')
MiDASGeno4[1370,] = c('HLA_A_0301', 'HLA_A_0301', 'HLA_C_0102', 'HLA_C_0401', 'HLA_B_3501', 'HLA_B_5101', 'HLA_DRB1_0101', 'HLA_DRB1_0101', 'HLA_DQA1_0101', 'HLA_DQA1_0101', 'HLA_DQB1_0501', 'HLA_DQB1_0501', 'HLA_DPA1_0103', 'HLA_DPA1_0103', 'HLA_DPB1_0301', 'HLA_DPB1_0402')
MiDASGeno4[1423,] = c('HLA_A_0101', 'HLA_A_0101', 'HLA_C_0501', 'HLA_C_0701', 'HLA_B_0801', 'HLA_B_4402', 'HLA_DRB1_0301', 'HLA_DRB1_0402', 'HLA_DQA1_0301', 'HLA_DQA1_0501', 'HLA_DQB1_0201', 'HLA_DQB1_0201', 'HLA_DPA1_0103', 'HLA_DPA1_0103', 'HLA_DPB1_0201', 'HLA_DPB1_0301')
MiDASGeno4[54,16] = 'HLA_DPB1_0402'

#### Check all columns contain the correct HLA type

MiDASGeno4 %>%
  filter(!grepl('HLA_A_', A_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_A_', A_2))

MiDASGeno4 %>%
  filter(!grepl('HLA_B_', B_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_B_', B_2)) 

MiDASGeno4 %>%
  filter(!grepl('HLA_C_', C_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_C_', C_2))

MiDASGeno4 %>%
  filter(!grepl('HLA_DRB1_', DRB1_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_DRB1_', DRB1_2))

MiDASGeno4 %>%
  filter(!grepl('HLA_DQA1_', DQA1_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_DQA1_', DQA1_2))

MiDASGeno4 %>%
  filter(!grepl('HLA_DQB1_', DQB1_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_DQB1_', DQB1_2))

MiDASGeno4 %>%
  filter(!grepl('HLA_DPA1_', DPA1_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_DPA1_', DPA1_2))

MiDASGeno4 %>%
  filter(!grepl('HLA_DPB1_', DPB1_1))
MiDASGeno4 %>%
  filter(!grepl('HLA_DPB1_', DPB1_2))


##### All columns OK #####

###Order columns

MiDASGeno4 = MiDASGeno4 %>%
  select(A_1,	A_2,	B_1,	B_2,	C_1,	C_2,	DPA1_1,	DPA1_2,	DPB1_1,	DPB1_2,	DQA1_1,	DQA1_2,	DQB1_1,	DQB1_2,	DRB1_1,	DRB1_2)
  

##### Adjust nomenclature to match MiDAS
MiDASGeno4 = MiDASGeno4 %>%
  mutate(across(everything(), str_replace_all, 'HLA_', '')) %>%
  mutate(across(everything(), str_replace_all, '_', '*'))
  

### Create function to add character
fun_insert = function(x, pos, insert){
  gsub(paste0("^(.{", pos,"})(.*)$"),
       paste0("\\1",insert,"\\2"),
       x)
}
# Ensure dataframe is a tibble rather than a matrix to be confident below function works as expected
MiDASGeno4=tibble(MiDASGeno4)

# Add : to second from last charcter of string
MiDASGeno4 = mutate(MiDASGeno4, across(everything(), ~ fun_insert(.x,str_length(.x)-2,':')))

# Export data table
write_tsv(MiDASGeno4, "HLA_files/MiDASGeno.txt")


