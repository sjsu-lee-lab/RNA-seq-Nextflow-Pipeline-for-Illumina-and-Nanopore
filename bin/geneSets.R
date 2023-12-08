#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
    
#Create -c argument that the user can use to specify their gene sets of interest
option_list <- list(optparse::make_option("--collections", 
  type = "character", 
  metavar = "C1,C2,C3,C4,C4,C5,C6,C7,C8,H", 
  help = "List of collections separated by commas"
))

# Parse the argument
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
collections <- unlist(strsplit(args$collections, ","))

# Get gene sets and combine them in one list
gene_sets <- lapply(collections, function(col) msigdbr(species = "Homo sapiens", category = col))
    
#Define the directory where the gene sets will be saved and save the gene sets to an RDS file in the specified directory
directory <- "gene_sets_dir"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(gene_sets, file.path(directory, 'gene_sets.rds'))
 

#Build a table dataframe of all collections with category (gs_cat), pathway (gs_name),and a list of its respective genes (gene_symbols) to be used in the python script for heatmaps 
table_all <- data.frame(gs_cat = character(), gs_name = character(), gene_list = I(list()))
for (gene_set in gene_sets) {
    table <- gene_set %>%
        group_by(gs_name, gs_cat) %>%
        summarise(gene_list = list(unique(gene_symbol))) %>%
        ungroup() %>% select(gs_cat, gs_name, gene_list)
    table$gene_list <- sapply(table$gene_list, toString)
    table_all <- rbind(table_all, table)
}
write.csv(table_all, file = "gene_sets.csv", row.names = FALSE)