setwd("~/PycharmProjects/microbe_masst/trees/microbe_masst_tree")

library(tidyverse)
library(taxize)

data <- read_csv("microbe_masst_table.csv")

# Extract unique NCBI IDs
unique_ncbi <- data %>% distinct(Taxa_NCBI) %>% dplyr::filter(Taxa_NCBI != "QC") %>% dplyr::filter(Taxa_NCBI != "Blank")
unique_ncbi$Taxa_NCBI <- as.numeric(unique_ncbi$Taxa_NCBI)

# Classify nodes tree
Sys.setenv(ENTREZ_KEY = "enter you key here")

microbe_lineage <- classification(unique_ncbi$Taxa_NCBI, db = 'ncbi', batch_size = 4) # it is slow (NCBI) and can return HTTP error - if so re run it
microbe_lineage_df <- do.call(rbind, microbe_lineage) %>% distinct(id, .keep_all = TRUE)
microbe_lineage_df_filter <- microbe_lineage_df %>% dplyr::filter(rank %in% c("no rank", "superkingdom", "kingdom", "phylum",
                                                                              "class", "order", "family", "genus", "subgenus",
                                                                              "species", "subspecies", "strain", "varietas"))

write_csv(microbe_lineage_df_filter, "list_ncbi_microbemasst.csv")
# I have saved the NCBI IDs as a list to generate a newick tree with ETE3 --> move to jupyter notebook



# Extracted NCBI IDs from ETE3 tree
tree_ids <- read_csv("Extracted_nodes_microbe.csv")
tree_ids$NCBI <- as.character(tree_ids$NCBI)

microbemasst_ids <- read_csv("list_ncbi_microbemasst.csv")
microbemasst_ids$id <- as.character(microbemasst_ids$id)

matching <- microbemasst_ids %>% left_join(tree_ids, by = c("id" = "NCBI")) 
# there is a problem with cellular organisms - no rank - 131567 --> check final JSON to manual fix

extra_nodes <- tree_ids %>% anti_join(microbemasst_ids, by = c("NCBI" = "id")) %>% dplyr::select(NCBI)
# these nodes will have to be removed from the final JSON cause it will otherwise affect the tree visualization levels

#write_csv(extra_nodes, "nodes_to_remove.csv") #remeber to remove NA from csv before reading it in the jupyter notebook

# Check number of available files per NCBI ID
availability <- data %>% group_by(Taxa_NCBI) %>% summarise(Presence = n()) %>% arrange(desc(Presence))
