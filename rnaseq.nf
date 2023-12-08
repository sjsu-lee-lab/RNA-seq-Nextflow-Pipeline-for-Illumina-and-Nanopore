#!/usr/bin/env nextflow

//use Nextflow DSL 2
nextflow.enable.dsl = 2
      
def helpMessage() {
  log.info """
   Usage:
    The typical command for running the pipeline is as follows:
        nextflow run rnaseq.nf

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
         --help                   This usage statement
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

if (params.nanopore) {
    log.info """\
           N A N O P O R E  R N A S E Q  P I P E L I N E    
         =================================================
         Results: ${params.outdir}
         Single-end reads : ${params.SEreads}
         Genome: ${params.genome}
         Annotation: ${params.annotation}
         Nanopore mode: ${params.nanopore}
         """
         .stripIndent()
} else {
    log.info """\
           I L L U M I N A  R N A S E Q  P I P E L I N E    
         =================================================
         Results: ${params.outdir}
         Paired-end reads : ${params.PEreads}
         Genome: ${params.genome}
         Annotation: ${params.annotation}
         Nanopore mode: ${params.nanopore}
         """
         .stripIndent()
}




//Quality control on raw Nanopore data
process NANOPLOT_RAW {
    tag "on $read_id"
    publishDir "$params.outdir/Nanopore/Nanoplot/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    path("${read_id}")
    
    when:
    params.nanopore
    
    script:
    """
    mkdir -p ${read_id}
    NanoPlot --fastq ${reads} -o ${read_id}
    """
}


//Quality control
process FASTQC {
    clusterOptions "--nodes=1"
    tag "on $read_id"
    if (params.nanopore) {
        label 'multi_long'
        publishDir "$params.outdir/Nanopore/FastQC/", mode:'copy'}
    else {
        label 'multi_short'
        publishDir "$params.outdir/Illumina/FastQC/", mode:'copy'}
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    path("fastqc_${read_id}")
    
    script:
    """
    mkdir -p fastqc_${read_id}
    fastqc -t ${task.cpus} -o fastqc_${read_id} -q ${reads}
    """
}


//Generating MultiQC report for raw data
process MULTIQC_RAW {
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/MultiQC/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/MultiQC/", mode:'copy'}
    
    input:
    path('*')
    
    output:
    path('multiqc_raw.html')
    
    script:
    """
    multiqc -n multiqc_raw.html .
    """
}


//Adapter trimming and filtering for short reads (quality_score>20;min_length=20;low complexity filter;threading;automatically detect adapters;json report)
process FASTP {
    label 'multi_short'
    clusterOptions "--nodes=1"
    tag "on $read_id"
    publishDir "$params.outdir/Illumina/Fastp/", mode:'copy', pattern: '*trimmed.fq.gz'
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    tuple val(read_id), path("${read_id}_{1,2}_trimmed.fq.gz"), emit: trimmed
    path("${read_id}.fastp.json"), emit: json
    
    when:
    !params.nanopore
    
    script:
    """
    fastp -i ${reads[0]} -o ${read_id}_1_trimmed.fq.gz -I ${reads[1]} -O ${read_id}_2_trimmed.fq.gz -q ${params.short_minQ} -l ${params.short_minLen} -y --thread ${task.cpus} --detect_adapter_for_pe --json ${read_id}.fastp.json
    """

}


//Adapter trimming for Nanopore
process PORECHOP {
    label 'multi_long'
    tag "on $read_id"
    publishDir "$params.outdir/Nanopore/Porechop/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    
    output:
    tuple val(read_id), path("${read_id}_trimmed.fastq")
    
    when:
    params.nanopore
    
    script:
    """
    porechop -t ${task.cpus} --discard_middle -i ${reads} -o ${read_id}_trimmed.fastq
    """
}


//Adapter filtering for Nanopore
process CHOPPER {
    label 'multi_long'
    tag "on $read_id"
    publishDir "$params.outdir/Nanopore/CHOPPER/", mode:'copy'
    
    input:
    tuple val(read_id), path(trimmed)
    
    output:
    tuple val(read_id), path("${read_id}_chopped.fastq")
    
    when:
    params.nanopore
    
    script:
    """
    cat ${trimmed} | chopper -q ${params.long_minQ} -l ${params.long_minLen} > ${read_id}_chopped.fastq
    """
}


//Quality control on Nanopore filtered data
process FASTQC_FILTERED {
    clusterOptions "--nodes=1"
    tag "on $read_id"
    publishDir "$params.outdir/Nanopore/FastQC/", mode:'copy'
    
    input:
    tuple val(read_id), path(filtered)
    
    output:
    path("${read_id}_filtered")
    
    when:
    params.nanopore
    
    script:
    """
    mkdir -p ${read_id}_filtered
    fastqc -o ${read_id}_filtered -q ${filtered}
    """
}


process NANOPLOT_FILTERED {
    tag "on $read_id"
    publishDir "$params.outdir/Nanopore/Nanoplot/", mode:'copy'
    
    input:
    tuple val(read_id), path(filtered)
    
    output:
    path("${read_id}_filtered")
    
    when:
    params.nanopore
    
    script:
    """
    mkdir -p ${read_id}_filtered
    NanoPlot --fastq ${filtered} -o ${read_id}_filtered
    """
}


process MULTIQC_FILTERED {
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/MultiQC/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/MultiQC/", mode:'copy'}
    
    input:
    path('*')
    
    output:
    path('multiqc_filtered.html')
    
    script:
    """
    multiqc -n multiqc_filtered.html .
    """
}


//Building index for Minimap2
process MINIMAP2_INDEX {
    memory '20 GB'
    publishDir "$params.outdir/Nanopore/Minimap2/Index", mode:'copy'
    
    input:
    path genome
    
    output:
    path("${genome.baseName}.mmi")
    
    when:
    params.nanopore
    
    script:
    """
    minimap2 -ax splice -d ${genome.baseName}.mmi ${genome}
    """

}


//Mapping Nanopore reads to the genome and converting sam file to coordinate-sorted bam file
process MINIMAP2_MAP {
    label 'multi_long'
    tag "on $read_id"
    publishDir "$params.outdir/Nanopore/Minimap2/", mode:'copy'
    
    input:
    tuple val(read_id), path(reads)
    path index
    
    output:
    path("map_${read_id}")
    tuple val(read_id), path("map_${read_id}/${read_id}_sorted.bam"), emit: sorted
    tuple val(read_id), path("map_${read_id}/${read_id}.bam.bai"), emit: index
    path("map_${read_id}/${read_id}.txt"), emit: stat
    
    when:
    params.nanopore
        
    script:
    """
    mkdir -p map_${read_id} 
    paftools.js gff2bed ${params.annotation} > annotation.bed
    minimap2 -t {task.cpus} -ax splice ${index} --junc-bed annotation.bed ${reads} | samtools view -hbo ${read_id}.bam
    samtools sort -o map_${read_id}/${read_id}_sorted.bam ${read_id}.bam
    samtools index map_${read_id}/${read_id}_sorted.bam map_${read_id}/${read_id}.bam.bai
    samtools stats map_${read_id}/${read_id}_sorted.bam > map_${read_id}/${read_id}.txt
    """
}


//Indexing the genome with STAR
process STAR_INDEX {
    label 'high'
    clusterOptions "--nodes=1"
    publishDir "$params.outdir/Illumina/STAR/", mode:'copy'
    
    output:
    path("Indicies")
    
    when:
    !params.nanopore
    
    script:
    """
    mkdir -p Indicies
    STAR --runThreadN 14 --runMode genomeGenerate --genomeDir Indicies --genomeFastaFiles ${params.genome} --sjdbGTFfile ${params.annotation} --sjdbOverhang ${params.overhang} --limitGenomeGenerateRAM 60000000000
    """
}


//Mapping filtered Illumina reads to genome generating a coordinate-sorted bam file
process STAR_MAP {
    label 'high'
    tag "on $read_id"
    publishDir "$params.outdir/Illumina/STAR/", mode:'copy'
    
    input:
    path indices
    tuple val(read_id), path(reads)
    
    when:
    !params.nanopore
    
    output:
    path("map_${read_id}")
    tuple val(read_id), path("map_${read_id}/${read_id}_Aligned.sortedByCoord.out.bam"), emit: sorted_bam
    path("map_${read_id}/${read_id}_Log.final.out"), emit: log_final
    tuple val(read_id), path("map_${read_id}/${read_id}.bam.bai"), emit: index
    
    
    script:
    """
    mkdir -p map_${read_id} 
    STAR --runThreadN 28 --genomeDir ${indices} --readFilesIn ${reads[0]} ${reads[1]}  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix map_${read_id}/${read_id}_
    samtools index map_${read_id}/${read_id}_Aligned.sortedByCoord.out.bam map_${read_id}/${read_id}.bam.bai
    """
}


process MULTIQC_ALIGNED {
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/MultiQC/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/MultiQC/", mode:'copy'}
    
    input:
    path('*')
    
    output:
    path('multiqc_aligned.html')
    
    script:
    """
    multiqc -n multiqc_aligned.html .
    """
} 


//Dealing with duplicates: AlignmentSummaryMetrics, MarkDuplicates, AlignmentSummaryMetrics
process DUPLICATES {
    tag "on $read_id"
    if (params.nanopore) {
        label 'multi_long'
        publishDir "$params.outdir/Nanopore/Duplicates/", mode:'copy'}
    else {
        label 'high'
        publishDir "$params.outdir/Illumina/Duplicates/", mode:'copy'}
    
    input:
    tuple val(read_id), path(bam)
   
   output:
    path("${read_id}_dup/${read_id}_metrics_before"), emit: premetrics
    path("${read_id}_dup/${read_id}_dedup.bam"), emit: dedup_bam
    path("${read_id}_dup/${read_id}.marked_dup_metrics.txt"), emit: dedupmetrics
    path("${read_id}_dup/${read_id}_metrics_after"), emit: postmetrics
    
    script:
    """
    mkdir -p ${read_id}_dup/${read_id}_metrics_before ${read_id}_dup/${read_id}_metrics_after
    picard -Xmx10G CollectMultipleMetrics I=${bam} O=${read_id}_dup/${read_id}_metrics_before/${read_id}_metrics_before R=${params.genome} PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=RnaSeqMetrics REF_FLAT=${params.refFlat}
    picard -Xmx10G MarkDuplicates I=$bam O=${read_id}_dup/${read_id}_dedup.bam M=${read_id}_dup/${read_id}.marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true
    picard -Xmx10G CollectMultipleMetrics I=${read_id}_dup/${read_id}_dedup.bam O=${read_id}_dup/${read_id}_metrics_after/${read_id}_metrics_after R=${params.genome} PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=RnaSeqMetrics REF_FLAT=${params.refFlat}
    """
}


process MULTIQC_DUPS {
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/MultiQC/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/MultiQC/", mode:'copy'}
    
    input:
    path('*')
    
    output:
    path('multiqc_dups.html')
    
    script:
    """
    multiqc -n multiqc_dups.html .
    """
}


//Quantifying reads (-p reads are paired, -B only reads with both ends successfully mapped are counted, -C chimeric reads are not counted, -L long-read mode)
process FEATURECOUNTS {
    if (params.nanopore) {
        label 'multi_long'
        publishDir "$params.outdir/Nanopore/FeatureCounts/", mode:'copy'}
    else {
        label 'multi_short'
        publishDir "$params.outdir/Illumina/FeatureCounts/", mode:'copy'}
    
    input:
    path(bams)
    
    output:
    path("counts.txt"), emit: counts
    path("counts.txt.summary"), emit: logs
    
    script:
    """
    if $params.nanopore
    then
        featureCounts -L -C -T 32 -t exon -g gene_id -a ${params.annotation} -s 0 -o counts.txt ${bams}
    else
        featureCounts -p -C -B -T 32 -t exon -g gene_id -a ${params.annotation} -s 0 -o counts.txt ${bams}
    fi
    """
}


process MULTIQC_COUNTS {
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/MultiQC/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/MultiQC/", mode:'copy'}
    
    input:
    path('*')
    
    output:
    path('multiqc_counts.html')
    
    script:
    """
    multiqc -n multiqc_counts.html .
    """
}


process DESEQ2 {
    conda "envs/r.yaml"
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/DESeq2/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/DESeq2/", mode:'copy'}
    
    input:
    path(counts)
    
    output:
    path("normalized_counts.csv"), emit: norm_counts
    path("DESeq_results.csv"), emit: DE_genes
    path("DESeq_summary.txt")
    path("Downregulated.csv")
    path("Upregulated.csv")
    path("GSEA_metric.rnk"), emit: ranked
    path("GSEA_shrunk.rnk")
    path("Rplots.pdf")
    
    script:
    """
    DeSeq.R --file ${counts}
    """
}

process FGSEA_SETS {
    conda "envs/r.yaml"
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/FGSEA/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/FGSEA/", mode:'copy'}
    
    output:
    path("gene_sets_dir"), emit: r_object
    path("gene_sets.csv"), emit: gene_sets
    
    script:
    """
    geneSets.R --collections ${params.collections}
    """
}

process FGSEA {
    debug true //to show echo in terminal
    conda "envs/r.yaml"
    label 'high'
    cpus '4'
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/FGSEA/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/FGSEA/", mode:'copy'}
    
    input:
    path(DE_genes)
    path(ranks)
    path(gene_sets_dir)
    
    output:
    path("Rplots.pdf"), optional: true
    path("*_top_pathways.csv"), optional: true, emit: top_pathways
    
    script:
    """
    fgsea.R --file ${ranks} --gs ${gene_sets_dir} --nperm ${params.nperm}
    """
}


process HEATMAPS {
    debug true //to show echo in terminal
    conda "envs/py.yaml"
    label 'high'
    cpus '4'
    if (params.nanopore) {publishDir "$params.outdir/Nanopore/FGSEA/Heatmaps/", mode:'copy'}
    else {publishDir "$params.outdir/Illumina/FGSEA/Heatmaps/", mode:'copy'}
    
    input:
    path(counts)
    path(DE_genes)
    path(gene_sets)
    path(top_pathways)
    
    output:
    path("heatmaps.pdf"), optional: true
    
    script:
    """
    lines="\$(cat $DE_genes | wc -l)"
    if [[ \$lines -gt 1 ]]
    then
        echo "DE genes detected. Heatmaps generated."
        python $projectDir/bin/heatmaps.py --counts ${counts} --g ${DE_genes} --gs ${gene_sets} --pathways ${top_pathways}
    else
        echo "No DE genes detected. Heatmaps omitted."
    fi
    """
}


workflow {
    if (params.nanopore) {
        //Convert the input from path into tuple(read_id, path)
        reads_ch = Channel
            .fromPath(params.SEreads, checkIfExists: true)
            .map { tuple( it.baseName, it ) }
    } else {
        reads_ch = Channel.fromFilePairs(params.PEreads, checkIfExists: true)
    }
    NANOPLOT_RAW(reads_ch)
    FASTQC(reads_ch)
    MULTIQC_RAW(FASTQC.out.collect())
    FASTP(reads_ch)
    PORECHOP(reads_ch)
    CHOPPER(PORECHOP.out)
    NANOPLOT_FILTERED(CHOPPER.out)
    FASTQC_FILTERED(CHOPPER.out)
    if (params.nanopore) {
        MULTIQC_FILTERED(FASTQC_FILTERED.out.collect())
    } else {MULTIQC_FILTERED(FASTP.out.json.collect())}
    MINIMAP2_INDEX(params.genome)
    MINIMAP2_MAP(CHOPPER.out, MINIMAP2_INDEX.out)
    STAR_INDEX()
    STAR_MAP(STAR_INDEX.out, FASTP.out.trimmed)
    if (params.nanopore) {
        MULTIQC_ALIGNED(MINIMAP2_MAP.out.stat.collect())
        DUPLICATES(MINIMAP2_MAP.out.sorted)
    } else {
        MULTIQC_ALIGNED(STAR_MAP.out.log_final.collect())
        DUPLICATES(STAR_MAP.out.sorted_bam)
        }
    MULTIQC_DUPS(DUPLICATES.out.premetrics.mix(DUPLICATES.out.postmetrics).mix(DUPLICATES.out.dedupmetrics).collect())
    FEATURECOUNTS(DUPLICATES.out.dedup_bam.collect())
    MULTIQC_COUNTS(FEATURECOUNTS.out.logs)
    DESEQ2(FEATURECOUNTS.out.counts)
    FGSEA_SETS()
    FGSEA(DESEQ2.out.DE_genes, DESEQ2.out.ranked, FGSEA_SETS.out.r_object)
    HEATMAPS(DESEQ2.out.norm_counts, DESEQ2.out.DE_genes, FGSEA_SETS.out.gene_sets, FGSEA.out.top_pathways.collect())
}


workflow.onComplete { 
    log.info ( workflow.success ? "\nSuccessfully Completed!" : "Oops .. Something went wrong: ${workflow.errorMessage}" )
}