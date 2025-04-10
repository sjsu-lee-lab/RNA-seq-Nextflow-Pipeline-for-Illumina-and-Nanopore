# RNASeq Pipeline for Illumina and Nanopore Samples using Nextflow 

This RNASeq pipeline was developed to evaluate the effect of walnut extract on Sp1-related pathways in triple negative breast cancer (TNBC) in humans (MDA-MB-231 cell line). It was built to run on the San Jose State University (SJSU) High Performance Cluster (HPC).

## Quick Setup

This pipeline requires Nextflow and Conda to execute its steps. To install Nextflow, you can run [`setup.sh`](setup.sh), that includes the necessary commands to install Nextflow and its dependencies:

```
./setup.sh </absolute/path/to/directory>
```
Then, you can add the installation directory to your PATH environment variable by adding the following line to your ~/.bashrc file:

```
export PATH="$PATH:</absolute/path/to/directory>"
```

Another way to install Nextflow is through Conda:

```
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```

To install the latest Conda, you can run the following:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

To check if Conda is working properly:
```
conda update conda
```


## Running the pipeline

To run the pipeline on Illumina samples:

```
sbatch Run_rnaseq.sh
```

To run the pipeline on Nanopore samples:

```
sbatch Run_rnaseq.sh --nanopore
```

## Calling help

```
nextflow run rnaseq.nf --help
```

## Output

Illumina results are stored in Illumina directory under the parent output directory `RNASeq_results` that is set with the `outdir` parameter as shown below. A timeline report, an execution report, and a DAG visualization are also generated for each run of the pipeline.

```
Illumina
├── DESeq2
│   ├── DESeq_results.csv
│   ├── DESeq_summary.txt
│   ├── Dowregulated.csv
│   ├── GSEA_metric.rnk
│   ├── GSEA_shrunk.rnk
│   ├── normalized_counts.csv
│   ├── Rplots.pdf
│   └── Upregulated.csv
├── Duplicates
├── Fastp
├── FastQC
├── FeatureCounts
├── FGSEA
│   ├── gene_sets_dir
│   ├── Heatmaps
│   │   └── heatmaps.pdf
│   ├── *_top_pathways.csv
│   ├── gene_sets.csv
│   ├── GSEA_shrunk.rnk
│   └── Rplots.pdf
├── MultiQC
│   ├── multiqc_aligned.html
│   ├── multiqc_counts.html
│   ├── multiqc_dups.html
│   ├── multiqc_filtered.html
│   └── multiqc_raw.html
├── STAR
├── dag.html
├── report.html
└── timeline.html
```

Similarly, Nanopore results are stored in `RNASeq_results/Nanopore` with `DESeq2`, `FGSEA`, and `MultiQC` sub-directories having the same content:

```
Nanopore
├── CHOPPER
├── DESeq2
├── Duplicates
├── FastQC
├── FeatureCounts
├── FGSEA
├── Minimap2
├── MultiQC
├── Nanoplot
├── Porechop
├── dag.html
├── report.html
└── timeline.html
```

## Custom parameters

This pipeline has the following custom parameters that you can set from the command line:

```
Optional arguments:
         --nanopore               Boolean to run pipeline in long-read mode. Default is short-read mode [false]
         --PEreads                Full path to short paired-end read files with Rex (example:*_{1,2}.fq.gz)
         --SEreads                Full path to long single-end read files with Rex (example:*.fastq)
         --genome                 Reference genome (full path required)
         --annotation             Gene annotation file (full path required)
         --refFlat                Gene annotations in refFlat form (full path required)
         --outdir                 Parent directory to place intermediate and final output files
         --[short/long]_minQ      Minimum quality score for filtering short/long reads
         --[short/long]_minLen    Minimum read length for short/long filtering
         --overhang               STAR's sjdboverhang (should be max(read length) - 1)
         --nperm                  FGSEA's number of permutations
         --collections            MSigDB collections for pathway enrichment analysis (C1,C2,C3,C4,C5,C6,C7,C8,H)
         --help                   The help usage statement
```

## Resources
Since the input data are part of a private research, their file paths are left empty in [`nextflow.config`](nextflow.config). The same thing goes to the paths of the genome and annotation files. Make sure you declare them from the command line using their corresponding parameters when running the pipeline. You can also add the file paths to their rescpective parameters directly in [`nextflow.config`](nextflow.config) before running the pipeline.

As for the custom scripts, the R scripts for DESeq2 and FGSEA and the python script for the heatmaps are placed in the [`bin`](bin/) directory. Specifically, DESeq.R was written while taking into account the project's sample's names. You might want to consider modifying the part that parses the sample's names (lines 23-25, 32-38) to match your samples.

The Conda environments used in this pipeline are specified in their respective YAML files in the [`envs`](envs/) directory.
