#!/usr/bin/env snakemake

#original author: Shaghayegh Soudi
#contributors: NA


configfile:
    "config.yaml"

SAMPLES = glob_wildcards(config['data']+"/{sample}.fastq")
### READS= ["R1","R2"] ## try read this way?


rule _symlink_fastq


rule _run_fastQC_raw_fastq
    input:
        fastq_1="{sample}.R1_CTAGTAGC-CTAGTAGC_fastq.zip",
        fastq_2="{sample}.R2_CTAGTAGC-CTAGTAGC_fastq.zip"
    output:
        html="{sample}_fastqc.html",
        zipfastq_1="{sample}_R1_fastqc.zip",
        zipfastq_2="{sample}_R2_fastqc.zip",
    log:
        stderr ="logs/indexes/{genome}/fastqc_raw_stderr.log"
    conda:
        "envs/fastqc.yaml"    
    params:
        threads=1 
    shell:
        "fastqc -o {output} -t {params.threads} {input.fastq_1} {input.fastq_2} 2> {log.stderr}"


rule _run_Trim_Galore
    input:
        forward="{sample}.R1_CTAGTAGC-CTAGTAGC_fastq.zip",
        reverse="{sample}.R2_CTAGTAGC-CTAGTAGC_fastq.zip"
    output:
        forward_trimmed="{sample}_R1_trimmed.fq.gz",
        reverse_trimmed="{sample}_R2_trimmed.fq.gz",
        report="{sample}_trimming_report.txt"
    params:
        adapter_options="--paired",  # Add additional options here if needed
        threads=1 
    conda:
        "envs/trim_galore.yaml"  # Path to the Trim Galore! conda environment YAML file
    shell:
        "trim_galore {params.adapter_options} -o ./ --paired {input.forward} {input.reverse} -o ./ --fastqc_args \"--threads {params.threads}\" 2> {output.report}"


rule _run_BWA_and_sort
    input:
        reference="indexes/{genome}/{genome}.fa"
        fastqR1=str(rules._run_Trim_Galore.output.forward_trimmed), 
        fastqR2=str(rules._run_Trim_Galore.output.reverse_trimmed)
    output:
        bam="{OUTPUT_DIR}/${sample}.sorted.bam"
    params:
        threads=4    
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    shell:
        "bwa mem -t {params.threads} {input.reference} {input.fastqR1} {input.fastqR2} | samtools sort -o {output}


rule _mark_duplicated
    input:
        bam=str(rules._run_BWA_and_sort.output.bam)
    output:
        bam_marked="aligned_reads_marked.bam",
        metrics="duplication_metrics.txt"
    log:
        "logs/mark_duplicates.log"
    params:
        ""  # optional params string
    shell:
        "java -jar ${SCRIPT}/picard.jar MarkDuplicates \
        I=${INPUT_DIR}/${SAMPLE_ID}.sorted.bam \
        O=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1.bam \
        M=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT"


rule _add_read_group
    input:
        "indexes/{genome}/{genome}.fa"
        "out/_val_1.fq.gz"
        "out/_val_2.fq.gz"
    output:
        "{OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    log:
        bam="logs/indexes/{genome}/Bisulfite_Genome.log"
    conda:
        "envs/samtools.yaml"  # Path to the Samtools conda environment YAML file    
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


rule _run_manta:
    input:
        bam=str(rules._add_read_group.output.bam),
        config="manta_config.ini"
    output:
        vcfs="{sample}_manta.vcf.gz",
        sv_bedpe="{sample}_manta_sv.bedpe.gz"
    params:
        threads=1,  # Number of threads to use, adjust as needed
        script = /home/users/shsoudi/emoding/envs/manta/bin
    conda:
        "envs/manta.yaml"  # Path to the Manta conda environment YAML file
    shell:
        "configManta.py --bam {input.bam} --config {input.config} --callRegions {params.threads} --callMemMb {params.memory * 1000} && \
         runWorkflow.py -m local -j {params.threads} -g {params.memory} && \
         mv results/variants/diploidSV.vcf.gz {output.vcfs} && \
         mv results/variants/diploidSV.bedpe.gz {output.sv_bedpe}"


rule _run_delly:
    input:
        bam="{sample}.sorted.bam",
        reference="genome.fa"
    output:
        vcf="{sample}.delly.vcf"
    params:
        threads=1,  # Number of threads to use, adjust as needed
        sv_type="DEL"  # Structural variant type (e.g., DEL, DUP, INV, etc.)
    conda:
        "envs/delly.yaml"  # Path to the Delly conda environment YAML file
    shell:
        "delly call -t {params.sv_type} -o {output.vcf} -g {input.reference} -x {input.bam}"  


rule run_smoove:
    input:
        bam="{sample}.sorted.bam",
        reference="genome.fa"
    output:
        vcf="{sample}.smoove.vcf"
    params:
        threads=1  # Number of threads to use, adjust as needed
    conda:
        "envs/smoove.yaml"  # Path to the Smoove conda environment YAML file
    shell:
        "smoove call -p {params.threads} -x --name {sample} -f {input.reference} -v {output.vcf} {input.bam}"               


rule run_gridss:
    input:
        bam="{sample}.sorted.bam",
        reference="genome.fa",
        gridss_config="gridss_config.ini"
    output:
        vcf="{sample}.gridss.vcf"
    params:
        threads=1  # Number of threads to use, adjust as needed
    conda:
        "envs/gridss.yaml"  # Path to the Gridss conda environment YAML file
    shell:
        "gridss --somatic -o {output.vcf} -r {input.reference} -b {input.bam} -v {input.gridss_config} --threads {params.threads}"


rule merge_sv_calls:
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