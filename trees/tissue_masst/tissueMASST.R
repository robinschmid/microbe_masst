setwd("~/OneDrive - University of California, San Diego Health/Projects/tissueMASST")

library(tidyverse)
library(jsonlite)
library(taxize)
library(ggpubr)

# Read latest ReDU file (downloaded on 1 Dec 2024)
redu <- read.delim("ReDU_Dec24.tsv") #687,585 files

# Keep only samples with MSMS data obtained from humans, mice, or rats
model_datasets <- redu %>% 
  dplyr::filter(str_detect(pattern = "Homo|Mus|Rat", NCBITaxonomy)) %>%
  dplyr::filter(str_detect(pattern = "10088|10090|10095|10114|10116|1383439|63221|9606", NCBITaxonomy)) %>%
  distinct(ATTRIBUTE_DatasetAccession) %>% 
  arrange(desc(ATTRIBUTE_DatasetAccession)) #1,667 datasets

redu_filter <- redu %>% 
  dplyr::filter(MS2spectra_count > 0) %>% #138,767 files
  dplyr::filter(ATTRIBUTE_DatasetAccession %in% model_datasets$ATTRIBUTE_DatasetAccession) %>% #48,615 files
  dplyr::filter(SampleTypeSub1 %in% c("biofluid", "blank_analysis", "blank_extraction", "blank_QC", 
                                      "missing value", "tissue")) %>% #47,595 files
  dplyr::filter(SampleType %in% c("animal", "blank_analysis", "blank_extraction", "blank_QC", "missing value")) %>% #47,585 files
  dplyr::filter(!(str_detect(pattern = "gas", ChromatographyAndPhase))) %>%
  dplyr::filter(str_detect(pattern = "mzML|mzXML", filename)) %>%
  dplyr::filter(str_detect(pattern = "10088|10090|10095|10114|10116|1383439|63221|9606|missing", NCBITaxonomy)) #45,065 files

tissuemasst <- redu_filter %>% dplyr::select(filename, ATTRIBUTE_DatasetAccession, SampleType, SampleTypeSub1,
                                             UBERONBodyPartName, HealthStatus, DOIDCommonName, NCBITaxonomy, 
                                             AgeInYears, LifeStage, BiologicalSex, UniqueSubjectID, USI, DataSource) %>%
  dplyr::mutate(UBERONBodyPartName = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                               str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                               str_detect(pattern = "skin", UBERONBodyPartName) ~ "skin",
                                               str_detect(pattern = "missing", UBERONBodyPartName) ~ SampleTypeSub1,
                                               TRUE ~ UBERONBodyPartName)) %>%
  dplyr::mutate(UBERONBodyPartName = gsub("missing value", "unknown", UBERONBodyPartName)) %>%
  dplyr::mutate(DOIDCommonName = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                           str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                           str_detect(pattern = "missing", DOIDCommonName) ~ HealthStatus,
                                           TRUE ~ DOIDCommonName)) %>%
  dplyr::mutate(DOIDCommonName = gsub("missing value", "unknown", DOIDCommonName)) %>%
  dplyr::mutate(NCBITaxonomy = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                         str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                         str_detect(pattern = "9606", NCBITaxonomy) ~ "Homo sapiens",
                                         str_detect(pattern = "Mus", NCBITaxonomy) ~ "Mus musculus",
                                         str_detect(pattern = "Rat", NCBITaxonomy) ~ "Rattus norvegicus",
                                         str_detect(pattern = "63221", NCBITaxonomy) ~ "Homo sapiens neanderthalensis",
                                         str_detect(pattern = "missing", NCBITaxonomy) ~ "unknown",
                                         TRUE ~ NCBITaxonomy)) %>%
  dplyr::filter(NCBITaxonomy != "unknown") %>%
  dplyr::mutate(BiologicalSex = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                          str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                          str_detect(pattern = "missing", BiologicalSex) ~ "unknown",
                                          TRUE ~ BiologicalSex)) %>%
  dplyr::filter(BiologicalSex != "asexual") %>%
  dplyr::mutate(LifeStage = case_when(str_detect(pattern = "blank_", SampleType) ~ "blank_qc",
                                      str_detect(pattern = "blank_", SampleTypeSub1) ~ "blank_qc",
                                      str_detect(pattern = "missing", LifeStage) ~ "unknown",
                                      str_detect(pattern = "^$", LifeStage) ~ "unknown",
                                      TRUE ~ LifeStage)) %>%
  dplyr::mutate(ID = paste(NCBITaxonomy, UBERONBodyPartName, DOIDCommonName, BiologicalSex, LifeStage, sep = "_"))

check <- tissuemasst %>% group_by(DOIDCommonName) %>% summarise(count = n()) %>% arrange(desc(count))

tissuemasst <- tissuemasst %>%
  dplyr::mutate(Filename = paste(ATTRIBUTE_DatasetAccession, sub(".*/", "", filename), sep = "/"))

#write_csv(tissuemasst, "tissue_masst_table.csv")

leaf_count <- tissuemasst %>% 
  group_by(NCBITaxonomy, UBERONBodyPartName, DOIDCommonName, BiologicalSex, LifeStage) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>%
  dplyr::mutate(ID = paste(NCBITaxonomy, UBERONBodyPartName, DOIDCommonName, BiologicalSex, LifeStage, sep = "_"))

#write.csv(leaf_count, "tissue_masst_id_count.csv")

