# Structural variant calling workflow with short and long read sequencing

This repository contains the pipelines for processing the illumina short read sequencing and pacbio and nanopore long read to detect structural variants and custom scripts for processing the vcf file and further downstream analysis and finally visualaizing the results. 
Some example plots are availabl in the plot direcory.
The variant calling workflow starts with quality control and alignment, similar to the other NGS applications. Alignments is folowed by alignment clean-up to prepapre data for variant calling.

There are three Sarcoma cell-lines files as example dataset sequenced by both long and short read sequencing.

# Detecting structural varisnts using Illumina short read sequencing


## Quality Control Steps
After demultiplexing the sequencing data, the workflow starts with quality control to ensure the sequences were high quality. two paired-end FASTQ files through a quality diagnostic and control pipeline.

The pipeline is 
1. Create base quality dianostic graphs.
2. Check reads for adaptor sequences.
3. Trim adaptor sequences
4. Trim poor quality bases.

Example plots are provided in the plots directory.
```Trim Galore``` is used for quality assesment and trimming the fastq files.
Here you can find more information about [Trim Galore](https://github.com/FelixKrueger/TrimGalore) 

Other recommended trimming programs:
- Trimmomatic
- Scythe

For Paired-end alignment pipeline uses ```BWA-mem```

```BWA-mem``` is recommended for high-quality queries as it is faster and more accurate and better performance.



Rscripts were also developed to process and visulaize vcf files. Rscripts are available in the scripts directory with example Rplots in the plots directory. 
 

# Detecting structural varisnts using PacBio and Oxford Nanopore long- read sequencing
