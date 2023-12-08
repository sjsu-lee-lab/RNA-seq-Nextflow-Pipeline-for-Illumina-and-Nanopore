#!/usr/bin/env python

import argparse
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

sns.set(color_codes=True)
sns.set(font_scale=0.8)
palette = sns.color_palette()

# Parse input
parser = argparse.ArgumentParser(description='Generate heatmaps of DE genes in enriched pathways')

# Add arguments with parameter names
parser.add_argument('--counts', type=str, help='Path to normalized counts table')
parser.add_argument('--g', type=str, help='Path to DESeq_results.csv')
parser.add_argument('--gs', type=str, help='Path to gene sets stored in gene_sets.csv')
parser.add_argument('--pathways', nargs="*", help='list of top pathways found in gene sets (*_top_pathways.csv)')

# Parse arguments
args = parser.parse_args()
counts = args.counts
genes = args.g
top_pathways = args.pathways
#print(top_pathways)

sets_df = pd.read_csv(args.gs)

#Split dataframe based on gs_cat and change its format to a list of dictionaries (collections) with pathway name as keys and list of genes as values 
sets_by_cat_df = sets_df.groupby("gs_cat")
gene_sets = []
for name, data in sets_by_cat_df:
    data['gene_list'] = data['gene_list'].apply(lambda x: x.split(', '))  # split at ", " to get rid of leading spaces
    gene_set = data.set_index('gs_name')['gene_list'].to_dict()
    gene_sets.append(gene_set)


# Get enriched pathways found by fgsea for each gene set to generate heatmaps for them only.

  # Extract pathway names from *_top_pathways.csv and store them in a list per collection (enriched_{category})
  # All lists are stored in one list "enriched"
categories = [top.rstrip("_top_pathways.csv") for top in top_pathways]
enriched = []
for item in range(len(top_pathways)):
    globals()[f'enriched_{categories[item]}'] = []
    globals()[f'enriched_{categories[item]}'] = pd.read_csv(top_pathways[item])
    globals()[f'enriched_{categories[item]}'] = globals()[f'enriched_{categories[item]}']["pathway"].tolist()
    enriched.append(globals()[f'enriched_{categories[item]}'])
#print(enriched_H)

  # Subset the gene sets dictionaries to keep only top pathways and their respective genes
enriched_pathways = []
for j in range(len(gene_sets)):
    new_dict = {}
    for path in gene_sets[j].keys():
        if path in enriched[j]:
            new_dict[path] = gene_sets[j].get(path)
    enriched_pathways.append(new_dict)
#print(len(enriched_pathways))
#print(enriched_pathways[0])


# Load datasets
data = pd.read_csv(counts)
#print(data.head())
dataset = data.set_index("Gene")
#print(len(dataset), dataset.head())


# Keep only DEGs in counts table
DE = pd.read_csv(genes)
DE = DE.set_index("Unnamed: 0")
#print(DE.head())
DEG = []
for gene in dataset.index:
    if gene in DE.index:
        DEG.append(gene)

DE_dataset = dataset.loc[DEG]
#print(len(DE_dataset), DE_dataset.head())


# Perform z-scale of normalized counts
z_df = DE_dataset.apply(lambda x: scipy.stats.zscore(x), axis=1)
#print(z_df.head())


# Color code samples by group
color_dict = {}

for sample in DE_dataset.columns:
    if "DMSO" in sample:
        color_dict[sample] = palette[1]
    else:
        color_dict[sample] = palette[4]

group_color = pd.Series(color_dict)
pdf_pages = PdfPages('heatmaps.pdf')

# Keep genes that are part of hallmark pathways only
for gene_set in enriched_pathways:
    # Remove genes that are not found in counts file from pathway gene lists to prevent thrown errors
    for pathway in gene_set.keys():
        not_found = []
        genes = gene_set.get(pathway)
        for gene in genes:
            if gene not in DE_dataset.index:
                not_found.append(gene)
        genes_in_counts = list(set(genes).difference(not_found))
        # Subset dataset based on genes involved in each pathway
        subset = z_df.loc[genes_in_counts]
        if len(subset) > 1:
            plot = sns.clustermap(subset, figsize=(8,9), cmap='RdBu_r', dendrogram_ratio=.05, col_colors=group_color)
            title = pathway + "\n Illumina Normalized Counts"
            plot.fig.suptitle(title, x=0.8, fontsize = "small")
            plot.ax_row_dendrogram.set_visible(False)
            plot.fig.subplots_adjust(right=0.7)
            plot.ax_cbar.set_position((0.85, .2, .03, .4))
            #plt.yticks(fontsize=9)
            pdf_pages.savefig(plot.fig)
            plt.close()
pdf_pages.close()