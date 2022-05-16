---
  title: "MiDAS quick start"

---

library("tidyverse")
library("kableExtra")
library("knitr")
library("midasHLA")
BiocManager::install("midasHLA")


## Reading phenotype  data
getwd()
#import phenotype data
Brainpheno= read_csv("Pheno_data/BrainPathology_tidyApoe.csv")
# Import HLA typing file as filepath
hla_calls_file <- "HLA_files/MiDASGeno.txt"

## Change HLA typing resolution as specified
#An error occurred because one HLA type had a space at end. To fix manually change
hla_calls <- readHlaCalls(hla_calls_file, resolution = 4)

## Compare allele frequencies,
# can change which HLA allele and population
freq_HLA <- getHlaFrequencies(hla_calls = hla_calls, compare = TRUE) %>%
  filter(Freq > 0.01)



freq_HLA_long <- tidyr::gather(
  data = freq_HLA,
  key = "population",
  value = "freq",
  "Freq",
  "USA NMDP African American pop 2",
  "USA NMDP Chinese",
  "USA NMDP European Caucasian",
  "USA NMDP Hispanic South or Central American",
  "USA NMDP Japanese",
  "USA NMDP North American Amerindian",
  "USA NMDP South Asian Indian",
  factor_key = TRUE
) %>% 
  filter(grepl("^C", allele))

plot_HLAfreq <-
  ggplot2::ggplot(data = freq_HLA_long, ggplot2::aes(x = allele, y = freq, fill = population)) +
  ggplot2::geom_bar(
    stat = "identity",
    position = ggplot2::position_dodge(0.7),
    width = 0.7,
    colour = "black"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = formattable::percent)

plot_HLAfreq

### HLA association analysis
##modify BrainPheno data to be compatible to analysis
Brainpheno2 = Brainpheno %>%
  add_column(ID = Brainpheno$IID) %>%
  select(ID,BraakStage) %>%
  mutate(across(everything(), ~replace(.x, .x == -9, NA)))

## Prepare MiDAS
HLA <- prepareMiDAS(
  hla_calls = hla_calls,
  colData = Brainpheno2,
  experiment = "hla_alleles"
)

# See alleles dataframe

df_alleles = as.data.frame(HLA)
view(df_alleles)


HLA <- HWETest(
  object = HLA,
  experiment = "hla_alleles",
  HWE_cutoff = 0.05 / 447,
  as.MiDAS = TRUE
)

HLA_model <- glm(BraakStage ~ term, data = HLA, family = gaussian())
HLA_results <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_alleles", 
  inheritance_model = "additive",
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

kableResults(HLA_results)

HLA_results_cond <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_alleles", 
  inheritance_model = "additive", 
  conditional = TRUE,
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

kableResults(HLA_results_cond, scroll_box_height = "200px")


##allelic divergence

HLA_het <- prepareMiDAS(
  hla_calls = hla_calls, 
  colData = Brainpheno2, 
  experiment = c("hla_het","hla_divergence")
)

# see as dataframe
df_HLA_het = as.data.frame(HLA_het)
view(df_HLA_het)

HLA_het_model <- glm(BraakStage ~ term, data=HLA_het, family=gaussian())

HLA_het_results <- runMiDAS(HLA_het_model, 
                            experiment = "hla_het", 
                            exponentiate = TRUE
)

kableResults(HLA_het_results)

HLA_div_results <- runMiDAS(HLA_het_model, 
                            experiment = "hla_divergence", 
                            exponentiate = TRUE
)

kableResults(HLA_div_results, scroll_box_height = "250px")

#### super types
Super <- prepareMiDAS(
  hla_calls = hla_calls,
  colData = Brainpheno2,
  experiment = "hla_supertypes"
)

# see as dataframe
df_Super = as.data.frame(Super)
view(df_Super)

# Correct for hardyWeinberg
Super <- HWETest(
  object = Super,
  experiment = "hla_supertypes",
  HWE_cutoff = 0.05 / 447,
  as.MiDAS = TRUE
)

#Model
HLA_model <- glm(BraakStage ~ term, data = Super, family = "gaussian")
HLA_results <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_supertypes", 
  inheritance_model = "additive",
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

#see results
kableResults(HLA_results)

HLA_results_cond <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_supertypes", 
  inheritance_model = "additive", 
  conditional = TRUE,
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

kableResults(HLA_results_cond, scroll_box_height = "200px")

#### Peptide binding sequence - Does not work with this data
Peptide <- prepareMiDAS(
  hla_calls = hla_calls,
  colData = Brainpheno2,
  experiment = "hla_g_groups"
)

# see as dataframe
df_Peptide = as.data.frame(Peptide)
view(df_Peptide)

Peptide <- HWETest(
  object = Peptide,
  experiment = "hla_g_groups",
  HWE_cutoff = 0.05 / 447,
  as.MiDAS = TRUE
)

HLA_model <- glm(BraakStage ~ term, data = Peptide, family = gaussian())
HLA_results <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_g_groups", 
  inheritance_model = "additive",
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

kableResults(HLA_results)

HLA_results_cond <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_supertypes", 
  inheritance_model = "additive", 
  conditional = TRUE,
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

kableResults(HLA_results_cond, scroll_box_height = "200px")

#### NK_ligand status
Ligand <- prepareMiDAS(
  hla_calls = hla_calls,
  colData = Brainpheno2,
  experiment = "hla_NK_ligands"
)

# see as dataframe
df_Ligand = as.data.frame(Ligand)
view(df_Ligand)

# Correct for hardyWeinberg
Ligand <- HWETest(
  object = Ligand,
  experiment = "hla_NK_ligands",
  HWE_cutoff = 0.05 / 447,
  as.MiDAS = TRUE
)

df_Ligand = as.data.frame(Ligand)
view(df_Ligand)

#Model
HLA_model <- glm(BraakStage ~ term, data = Ligand, family = gaussian())
HLA_results <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_NK_ligands", 
  inheritance_model = "additive",
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

#see results
kableResults(HLA_results)

HLA_results_cond <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_NK_ligands", 
  inheritance_model = "additive", 
  conditional = TRUE,
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "fdr", 
  exponentiate = TRUE
)

kableResults(HLA_results_cond, scroll_box_height = "200px")
