#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

option_list = list(optparse::make_option("--file",
                                        type = 'character',
                                        help = "path to count file"))

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

#Read counts file with "Geneid" column as the row names
counts_table <- read.table(opt$file, header = TRUE, sep = "\t", row.names = 1)

#Remove first 5 columns (chr,start,end,strand,length)
counts_table <- counts_table[,-c(1:5)]

#Remove _dedup and _231 from sample names
colnames(counts_table) <- sub("_dedup.bam", "", colnames(counts_table))
if(any(grepl("_231", colnames(counts_table)))){
    colnames(counts_table) <- sub("_231", "", colnames(counts_table))}

#Drop unwanted samples
#if(any(grepl("WE_A", colnames(counts_table)))){
#counts_table = subset(counts_table, select = -WE_A )}

#Create sample information
samples <- colnames(counts_table)
if(any(grepl("_", samples))){
    sample_group <- strsplit(samples, "_") #to split condition from sample name
    group <- rapply(sample_group, function(x) head(x, 1)) #to get condition (1st element)
} else {
    group <- substring(samples, 3) #to get rid of first 2 characters
}
colData <- data.frame(condition=factor(group))
                

#Build dataset, run differential expression, and get results with an adjusted p-value cut-off of 0.05
dds <- DESeqDataSetFromMatrix(countData = counts_table, colData = colData, design = ~ condition)
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "WE", "DMSO"), alpha = 0.05)

#Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- cbind(rownames(normalized_counts), normalized_counts) #to make the index the first column of dataframe
colnames(normalized_counts)[1] <- "Gene" #to rename the first column
                
#Filter DE genes
DE_genes <- as.data.frame(results) %>%
    filter(abs(log2FoldChange) >= 0, padj<=0.05) %>% 
        arrange(padj)
                
#Upregulated genes
up_genes <- DE_genes %>%
  filter(log2FoldChange > 0) %>%
    arrange(padj, desc(log2FoldChange))

#Downregulated genes
down_genes <- DE_genes %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj, log2FoldChange)

#Prepare ranked list of genes sorted by the metric formula "-log10({p value}) * sign({Fold Change})" for gsea
df <- as.data.frame(results)
df <- df[!is.na(df$pvalue),] #get rid of NA in p-value
df <- df[!is.na(df$padj),] #get rid of NA in padj
df$pvalue[df$pvalue == 0] <- 1e-50  #change p-value of 0 to a very low number to prevent Inf in metrics
df$signFC <- sign(df$log2FoldChange)
df$logP <- -log10(df$pvalue)
df$metric <- df$logP * df$signFC
ranked_metric <- data.frame(row.names(df), df$metric)
colnames(ranked_metric) = c("Gene_ID", "metric")
ranked_metric <- ranked_metric[order(-ranked_metric$metric), ]
                
#Shrink LFC and prepare ranked list of shrunken LFC values for gsea
res_shrunk <- lfcShrink(dds, coef=resultsNames(dds)[2], type="ashr") #coef should be in resultsNames(dds)
df_shrunk <- as.data.frame(res_shrunk)
df_shrunk <- df_shrunk[!is.na(df_shrunk$pvalue),] #get rid of NA in p-value
df_shrunk <- df_shrunk[!is.na(df_shrunk$padj),] #get rid of NA in padj
shrunk_gsea <- data.frame(rownames(df_shrunk), df_shrunk$log2FoldChange)
colnames(shrunk_gsea) = c("Gene_ID", "log2FoldChange")
shrunk_gsea <- shrunk_gsea[order(-shrunk_gsea$log2FoldChange), ]

                
#Write results to output file
sink(file="DESeq_summary.txt")
summary(results)
sink()
write.table(normalized_counts, file="normalized_counts.csv", sep=",", quote=F, row.names=FALSE)
write.table(DE_genes, file="DESeq_results.csv", sep=",", quote=F, col.names=NA)
write.table(up_genes, file="Upregulated.csv", sep=",", quote=F, col.names=NA)
write.table(down_genes, file="Downregulated.csv", sep=",", quote=F, col.names=NA)
write.table(ranked_metric, file="GSEA_metric.rnk", sep="\t", quote=F, row.names=FALSE)
write.table(shrunk_gsea, file="GSEA_shrunk.rnk", sep="\t", quote=F, row.names=FALSE)              

                
#MA plot of normal results
plotMA(dds, ylim=c(-6,6), main="MA plot with no shrinkage")
                
#MA plot of shrinked LFC reducing noise
plotMA(res_shrunk, ylim=c(-6,6), main="MA plot with shrinkage")

                
#Transform dataset
rld <- rlog(dds, blind = TRUE) #blind=true: transformation unbiased to sample condition information

#PCA plot
p <- plotPCA(rld, intgroup="condition") #intgroup should specify columns of colData(dds)
print(p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = factor(colnames(dds))), size = 3)  + geom_point()) #requires ggplot2 and ggrepel

#Sample-sample comparison heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color=heat.colors, main = "Sample-to-Sample correlation heatmap",fontsize = 10, fontsize_row = 10, height = 20)

#Dispersion estimates
plotDispEsts(dds, main = "Dispersion estimates") #on the DESeq dataset           
                
#P-value distribution
hist(results$pvalue,main="Pvalue distribution",xlab="P-values")
                    
#P-adjusted distribution
hist(results$padj,main="P-adjusted distribution",xlab="P-adjusted")