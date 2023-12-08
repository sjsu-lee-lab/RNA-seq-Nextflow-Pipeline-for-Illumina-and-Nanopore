#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(msigdbr))
#suppressPackageStartupMessages(library(tibble)) #to use for changing index to column in dataframes (no longer needed)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))

pdf.options(width=18 , height=8)

fgsea_results <- function(ranks, gene_set, set_names, nperm){
	#Group gene IDs by gene set names
	gene_list <- split(x = gene_set$gene_symbol, f = gene_set$gs_name)
	#Run fgsea
	res <- fgsea(gene_list, ranks, minSize=15, maxSize=500, nPermSimple = nperm)
	#Plot enrichment table
	topUp <- res %>% filter(ES > 0) %>% arrange(padj) %>% head(n=10)
	topDown <- res %>% filter(ES < 0) %>% arrange(padj) %>% head(n=10)
	top <- bind_rows(topUp, topDown) %>% arrange(-NES)
    ttl <- paste(as.character(gene_set$gs_cat[1]), "_top_pathways.csv", sep = "")
	top_paths <- top[,-c(2:8)]
	write.csv(top_paths,file=ttl,row.names=F)
	title <- paste("Enrichment plot for the", as.character(set_names[gene_set$gs_cat[1]]), "gene set")
	#plot.new() not needed with grid.arrange
	grid.arrange(left = title, plotGseaTable(gene_list[top$pathway], ranks, res, gseaParam=0.5, colwidths = c(5,3,0.8,0.8,0.8)))
}

#Get DEG ranked list and list of gene sets
option_list = list(optparse::make_option("--file",
                                        type = 'character',
                                        help = "path to ranked list file"))
option_list = c(option_list, list(optparse::make_option("--gs",
                                        type = 'character',
                                        help = "directory holding list of gene sets as R object")))
option_list = c(option_list, list(optparse::make_option("--nperm",
                                        type = 'integer',
                                        help = "Number of permutations for FGSEA")))

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

ranks <- read.table(opt$file, header = TRUE, colClasses = c("character", "numeric"))
#assign gene IDs to metrics
ranks <- setNames(ranks$metric, ranks$Gene_ID)

gene_sets <- readRDS(file.path(opt$gs, 'gene_sets.rds'))

nperm = opt$nperm

set_names <- c(
	"H"="Hallmark (H)",
	"C2"="Curated (C2)",
	"C3"="Regulatory (C3)",
	"C4"="Computational (C4)",
	"C5"="Ontology (C5)",
	"C6"="Oncogenic (C6)",
	"C7"="Immunologic (C7)"
)

# Run analysis on each of the gene sets
for (gene_set in gene_sets) {
	fgsea_results(ranks, gene_set, set_names, nperm)
}


#To save the results in a text format data:table::fwrite function can be used:
#fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))