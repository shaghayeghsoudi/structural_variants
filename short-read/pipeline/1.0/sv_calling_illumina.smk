#!/usr/bin/env snakemake

#original author: Shaghayegh Soudi
#contributors: NA


configfile:
    "config.yaml"

SAMPLES = glob_wildcards(config['data']+"/{sample}.fastq")




rule run_fastQC_raw_fastq
    input:
        "indexes/{genome}/{genome}.fa"
    output:
        directory("indexes/{genome}/Bisulfite_Genome")
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        "v2.6.0/bio/bismark/bismark_genome_preparation"



rule run_Trim_Galore
    input:
        "indexes/{genome}/{genome}.fa"
    output:
        directory("indexes/{genome}/Bisulfite_Genome")
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        "v2.6.0/bio/bismark/bismark_genome_preparation"


rule run_BWA_and_sort
    input:
        "indexes/{genome}/{genome}.fa"
        "out/_val_1.fq.gz"
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        bwa mem -t 10 ${FASTA} ${SAMPLE_NAME_R1} ${SAMPLE_NAME_R2} | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam


rule mark_duplicated
input:
        "indexes/{genome}/{genome}.fa"
        "out/_val_1.fq.gz"
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        bwa mem -t 10 ${FASTA} ${SAMPLE_NAME_R1} ${SAMPLE_NAME_R2} | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam



rule add_read_group
input:
        "indexes/{genome}/{genome}.fa"
        "out/_val_1.fq.gz"
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        bwa mem -t 10 ${FASTA} ${SAMPLE_NAME_R1} ${SAMPLE_NAME_R2} | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam



rule samtools_index
input:
        "indexes/{genome}/{genome}.fa"
        "out/_val_1.fq.gz"
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        bwa mem -t 10 ${FASTA} ${SAMPLE_NAME_R1} ${SAMPLE_NAME_R2} | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam


rule manta
input:
        "indexes/{genome}/{genome}.fa"
        "out/_val_1.fq.gz"
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        bwa mem -t 10 ${FASTA} ${SAMPLE_NAME_R1} ${SAMPLE_NAME_R2} | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam
