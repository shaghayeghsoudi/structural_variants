#!/usr/bin/env snakemake

#original author: Shaghayegh Soudi
#contributors: NA


configfile:
    "config.yaml"

SAMPLES = glob_wildcards(config['data']+"/{sample}.fastq")




rule run_fastQC_raw_fastq
    input:
        fastq="{sample}.fastq"
    output:
        html="{sample}_fastqc.html",
        zip="{sample}_fastqc.zip"
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    conda:
        "envs/fastqc.yaml"    
    params:
        threads=1 
    shell:
        "fastqc -o {output} -t {params.threads} {input.fastq}"



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
        reference="indexes/{genome}/{genome}.fa"
        fastq1="sample.fastq"
        fastq2="sample.fastq"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    params:
        threads=4    
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    shell:
        "bwa mem -t {params.threads} {input.reference} {fastq1} {fastq2} | samtools sort -o {output}


rule mark_duplicated
input:
        bam="aligned_reads.bam"
    output:
        bam_marked="aligned_reads_marked.bam",
        metrics="duplication_metrics.txt"
    log:
        "logs/mark_duplicates.log"
    params:
        ""  # optional params string
   shell:
        "java -jar ${SCRIPT}/picard.jar AddOrReplaceReadGroups \
        I=${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam \
        O=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1.RG.bam \
        RGID=${SAMPLE_ID} \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=${SAMPLE_ID} \
        RGSM=${SAMPLE_ID} 2> ${SAMPLE_ID}.log"


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
