#!/usr/bin/env snakemake

#original author: Shaghayegh Soudi
#contributors: NA


configfile:
    "default.yaml"

DIRS,SAMPLES,READS,I7,I5 = glob_wildcards(config['data']+"{dir}/{sample}.{read}_{i7}-{i5}_fastq.zip")


rule _all:    ### Final TO DO, adjust rule all in the end based on what utputs are final
    input:
        expand("results/mapped/{dir}/{sample}.XX", zip, dir=DIRS, sample=SAMPLES, genome=GENOMES),
        expand("results/infer_sex/{dir}/{sample}.normal.sex.txt",zip, dir=DIRS, sample=SAMPLES),
        expand("results/infer_sex/{dir}/{sample}.normal.sex.txt",zip, dir=DIRS, sample=SAMPLES)


rule _symlink_fastq:
    input:
         fastq_1=config['data']+"{dir}/{sample}.R1_{i7}-{i5}_fastq.zip",
         fastq_2=config['data']+"{dir}/{sample}.R2_{i7}-{i5}_fastq.zip"
    output:
         fastq_1="00-input/{dir}/{sample}.R1_{i7}-{i5}_fastq.zip",
         fastq_2="00-input/{dir}/{sample}.R2_{i7}-{i5}_fastq.zip"
    shell: 
        "ln -s {input}.fastq_1 {output}.fastq_2 && ln -s {input}.fastq_2 {output}.fastq_2"     
               

rule _run_fastQC_raw_fastq
    input:
        fastq_1= str(rules._symlink_fastq.output.fastq_1),
        fastq_2= str(rules._symlink_fastq.output.fastq_2)
    output:
        html_1="{sample}_{i7}-{i5}_fastqc.html",
        html_2="{sample}_{i7}-{i5}_fastqc.html",
        zipfastq_1="{sample}_{i7}-{i5}_R1_fastqc.zip",
        zipfastq_2="{sample}_{i7}-{i5}_R2_fastqc.zip"
    log:
        stderr_fastq ="logs/indexes/{sample}_fastqc_raw_stderr.log"
    conda:
        "envs/fastqc.yaml"    
    params:
        threads=3
        out_dir = "/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV/fastQC_raw"
    shell:
        "fastqc -o {params.out_dir} -t {params.threads} {input.fastq_1} {input.fastq_2} 2> {log.stderr_fastq}"


rule _run_Trim_Galore
    input:
        forward="{sample}.R1_CTAGTAGC-CTAGTAGC_fastq.zip",
        reverse="{sample}.R2_CTAGTAGC-CTAGTAGC_fastq.zip"
    output:
        R1_trimmed="{sample}_R1_trimmed.fq.gz",
        R2_trimmed="{sample}_R2_trimmed.fq.gz",
        report="{sample}_trimming_report.txt"
    log:
        stderr ="logs/indexes/{sample}_fastqc_raw_stderr.log"    
    params:
        adapter_options="--paired",  
        threads=3,
        out_dir = "/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV/trimmed_fastq",
        cutadapt_script="/home/users/shsoudi/emoding/envs/trim-galore/bin/cutadapt"
    conda:
        "envs/trim_galore.yaml"  
    shell:
        """
        trim_galore --gzip \
        --path_to_cutadapt {params.cutadapt_script} \
        {params.adapter_options} \
        --output_dir {params.out_dir} \
        {input.forward} {input.reverse}  \
        --fastqc_args \
        --threads {params.threads} \
        2> {log.stderr}
        """


rule _run_BWA_and_sort
    input:
        reference="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa"
        fastqR1=str(rules._run_Trim_Galore.output.R1_trimmed), 
        fastqR2=str(rules._run_Trim_Galore.output.R2_trimmed)
    output:
        bam="/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV/alignment-BWA/{sample}/{sample}.sorted.bam"
    conda:
        "envs/bwa.yaml"    
    params:
        threads=10,
    log:
        stderr="logs/indexes/BWA_stderr.log"
    shell:
        """
        bwa mem -t {params.threads} \
        {input.reference} {input.fastqR1} {input.fastqR2} | samtools sort -o {output}.bam 2> {log.stderr}
        """


rule _mark_duplicated
    input:
        bam=str(rules._run_BWA_and_sort.output.bam)
    output:
        bam_marked="aligned_reads_marked.bam",
        metrics="duplication_metrics.txt"
    log:
        stderr="logs/mark_duplicates.log"
    params:
        script="/home/users/shsoudi/emoding/envs/picard/share/picard-2.27.5-0"
    shell:
        """
        java -jar {params.script}/picard.jar MarkDuplicates \
        I={input.bam} \
        O=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1.bam \
        M=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT 2> {log.stderr}
        """


rule _add_read_group
    input:
        "indexes/{genome}/{genome}.fa",
        "out/_val_1.fq.gz",
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        stderr="logs/indexes/{genome}/Bisulfite_Genome.log"
    conda:
        "envs/samtools.yaml"  # Path to the Samtools conda environment YAML file    
    params:
        ""  # optional params string
    shell:
        """
        java -jar ${SCRIPT}/picard.jar AddOrReplaceReadGroups \
        I=${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam \
        O=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1.RG.bam \
        RGID=${SAMPLE_ID} \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=${SAMPLE_ID} \
        RGSM=${SAMPLE_ID} 2> ${log.stderr}
        """


rule _samtools_index
    input:
        bam=str(rules._add_read_group.output.bam)
    output:
        bai="{sample}.sorted.bam.bai"
    params:
        threads=1  # Number of threads to use, adjust as needed
    conda:
        "envs/samtools.yaml"  # Path to the Samtools conda environment YAML file
    shell:
        "samtools index -@ {params.threads} {input.bam}"


rule configure_manta:
    input:
        bam=str(rules._add_read_group.output.bam),
        bai=str(rules._add_read_group.output.bai),
        reference="path/to/reference.fasta"
    output:
        config_dir="path/to/manta_config_dir"
    log:
        stderr="logs/indexes/{genome}/manta.log"     
    params:
        manta_config="path/to/manta_config.json"
    conda:
        "../envs/manta.yaml"  # Path to the Manta conda environment YAML file    
    shell:
        """
        mkdir -p {output.config_dir}
        manta/bin/configManta.py --bam {input.bam} --referenceFasta {input.reference} --runDir {output.config_dir} --config {params.manta_config}
        """


rule _run_manta:
    input:
        bam=str(rules._add_read_group.output.bam),
    output:
        vcfs="{sample}_manta.vcf.gz",
        sv_bedpe="{sample}_manta_sv.bedpe.gz"
    log:
        stderr="logs/indexes/{genome}/manta.log"    
    params:
        config.callers.manta.threads
        script = /home/users/shsoudi/emoding/envs/manta/bin    ### TO DO: adjust script in config
    conda:
        "../envs/manta.yaml"  # Path to the Manta conda environment YAML file
    shell:
        """
        configManta.py --bam {input.bam} --config {input.config} --callRegions {params.threads} --callMemMb {params.memory * 1000} && \
         runWorkflow.py -m local -j {params.threads} -g {params.memory} && \
         mv results/variants/diploidSV.vcf.gz {output.vcfs} && \
         mv results/variants/diploidSV.bedpe.gz {output.sv_bedpe} 2> {log.stderr}
         """


rule _run_delly:
    input:
        bam="str(rules._add_read_group.output.bam),
        reference="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa"
    output:
        vcf="{sample}.delly.vcf"
    log:
        stderr="logs/indexes/{genome}/manta.log"        
    params:
        config.callers.delly.threads
    conda:
        "envs/delly.yaml"  # Path to the Delly conda environment YAML file
    shell:
        """
        for name in DEL DUP INV TRA INS; 
        do delly call -t ${name} -o ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.bcf -g ${FASTA} ${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam  &&
        bcftools view ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.bcf > ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.vcf 2> {log.stderr} 
        done
        """ 

rule _delly_merge:


rule _run_smoove:
    input:
        bam="{sample}.sorted.bam",
        reference="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa"
    output:
        vcf="{sample}.smoove.vcf"
    params:
        threads=1  # Number of threads to use, adjust as needed
    conda:
        "envs/smoove.yaml"  # Path to the Smoove conda environment YAML file
    shell:
        "smoove call -p {params.threads} -x --name {sample} -f {input.reference} -v {output.vcf} {input.bam}"               


rule _run_gridss:
    input:
        bam="{sample}.sorted.bam",
        reference="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa",
        gridss_config="gridss_config.ini"
    output:
        vcf="{sample}.gridss.vcf"
    params:
        threads=1  # Number of threads to use, adjust as needed
    conda:
        "envs/gridss.yaml"  # Path to the Gridss conda environment YAML file
    shell:
        "gridss --somatic -o {output.vcf} -r {input.reference} -b {input.bam} -v {input.gridss_config} --threads {params.threads}"


rule _merge_sv_calls:
    input:
        vcf_files=expand("{sample}.vcf", sample=SAMPLES),
    output:
        merged_vcf="{output_dir}/merged_sv_calls.vcf"
    params:
        config="config.txt"  # Survivor configuration file
    conda:
        "envs/survivor.yaml"  # Path to the Survivor conda environment YAML file
    shell:
        "SURVIVOR merge {input.vcf_files} {params.config} 1000 1 1 1 1 {output.merged_vcf}"


rule visualization
rule run_r_script:
    input:
        "input.txt"
    output:
        "output.txt"
    script:
        "my_analysis.R"