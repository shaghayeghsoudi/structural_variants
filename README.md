# Detecting structural varisnts using Illumina short read sequencing
This repository contains the pipeline for processing the illumina short read sequencing to detect structural variants and scripts for processing the vcf file and downstream analysis and finally visualaizing the results. 
Some example plots are availabl in the plot direcory.

We sequenced three Sarcoma cell-lines (GCT, RD and SW), using Illumina paired-end sequencing.

## Quality Control Steps
After demultiplexing the sequencing data, the first step of analysis was to ensure the sequences were high quality. We ran each of the three lines. two paired-end FASTQ files through a quality diagnostic and control pipeline.

The pipeline is 
1. Create base quality dianostic graphs.
2. Check reads for adaptor sequences.
3. Trim adaptor sequences
4. Trim poor quality bases.

```Trim Galore``` is used for quality assesment and trimming the fastq files.
Here you can find more information about [Trim Galore](https://github.com/FelixKrueger/TrimGalore) 
