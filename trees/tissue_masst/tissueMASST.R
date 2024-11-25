setwd("~/OneDrive - University of California, San Diego Health/Projects/tissueMASST")

library(tidyverse)
library(jsonlite)
library(taxize)
library(ggpubr)

# Read latest ReDU file (downloaded on 25 Nov 2024)
redu <- read.delim("ReDU_Nov24.tsv") #681,434 files

# Keep only samples with MSMS data obtained from humans, mice, or rats
model_datasets <- redu %>% 
  dplyr::filter(str_detect(pattern = "Homo|Mus|Rat", NCBITaxonomy)) %>%
  distinct(ATTRIBUTE_DatasetAccession) %>% 
  arrange(desc(ATTRIBUTE_DatasetAccession)) #1,669 datasets

redu_filter <- redu %>% 
  dplyr::filter(MS2spectra_count > 0) %>% #135,258 files
  dplyr::filter(ATTRIBUTE_DatasetAccession %in% model_datasets$ATTRIBUTE_DatasetAccession) %>% #54,807 files
  dplyr::filter(SampleTypeSub1 %in% c("biofluid", "blank_analysis", "blank_extraction", "blank_QC", 
                                      "missing value", "tissue")) %>% #53,787 files
  dplyr::filter(SampleType %in% c("animal", "blank_analysis", "blank_extraction", "blank_QC", "missing value")) %>% #48,225 files
  dplyr::filter(str_detect(pattern = "9606|10088|10090|10116|10095|10114|63221|missing", NCBITaxonomy)) #45,876 files

tissuemasst <- redu_filter %>% dplyr::select(filename, ATTRIBUTE_DatasetAccession, SampleType, SampleTypeSub1,
                                             UBERONBodyPartName, HealthStatus, DOIDCommonName, NCBITaxonomy, 
                                             AgeInYears, LifeStage, BiologicalSex, UniqueSubjectID, USI, DataSource) %>%
  dplyr::mutate(UBERONBodyPartName = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                               str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                               str_detect(pattern = "skin", UBERONBodyPartName) ~ "skin",
                                               str_detect(pattern = "missing", UBERONBodyPartName) ~ SampleTypeSub1,
                                               TRUE ~ UBERONBodyPartName)) %>%
  dplyr::mutate(UBERONBodyPartName = gsub("missing value", "unknonw", UBERONBodyPartName)) %>%
  dplyr::mutate(DOIDCommonName = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                           str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                           str_detect(pattern = "missing", DOIDCommonName) ~ HealthStatus,
                                           TRUE ~ DOIDCommonName)) %>%
  dplyr::mutate(DOIDCommonName = gsub("missing value", "unknonw", DOIDCommonName)) %>%
  dplyr::mutate(NCBITaxonomy = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                         str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                         str_detect(pattern = "9606", NCBITaxonomy) ~ "Homo sapiens",
                                         str_detect(pattern = "Mus", NCBITaxonomy) ~ "Mus musculus",
                                         str_detect(pattern = "Rat", NCBITaxonomy) ~ "Rattus norvegicus",
                                         str_detect(pattern = "63221", NCBITaxonomy) ~ "Homo sapiens neanderthalensis",
                                         str_detect(pattern = "missing", NCBITaxonomy) ~ "unknonw",
                                         TRUE ~ NCBITaxonomy)) %>%
  dplyr::filter(NCBITaxonomy != "unknonw") %>%
  dplyr::mutate(BiologicalSex = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                          str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                          str_detect(pattern = "missing", BiologicalSex) ~ "unknonw",
                                          TRUE ~ BiologicalSex)) %>%
  dplyr::filter(BiologicalSex != "asexual") %>%
  dplyr::mutate(LifeStage = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                      str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                      str_detect(pattern = "missing", LifeStage) ~ "unknonw",
                                      str_detect(pattern = "^$", LifeStage) ~ "unknonw",
                                      TRUE ~ LifeStage)) %>%
  dplyr::mutate(ID = paste(UBERONBodyPartName, DOIDCommonName, NCBITaxonomy, BiologicalSex, LifeStage, sep = "_"))

check <- tissuemasst %>% group_by(DOIDCommonName) %>% summarise(count = n()) %>% arrange(desc(count))

write_csv(tissuemasst, "tissue_masst_table.csv")

leaf_count <- tissuemasst %>% 
  group_by(UBERONBodyPartName, DOIDCommonName, NCBITaxonomy, BiologicalSex, LifeStage) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>%
  dplyr::mutate(ID = paste(UBERONBodyPartName, DOIDCommonName, NCBITaxonomy, BiologicalSex, LifeStage, sep = "_"))

write.csv(leaf_count, "tissue_masst_id_count.csv")
