#!/bin/bash

### activate conda
source /home/shsoudi/anaconda3/etc/profile.d/conda.sh
source activate /home/shsoudi/anaconda3/envs/gridds

# Path to the reference genome
reference_genome=/home/shsoudi/oak-indices/hg19.fa

# Base directory containing BAM files (if sample file only contains sample names)
bam_dir="/home/shsoudi/oak/cell-lines/alignment-BWA/dedupped_BAMs_picard"
normal_bam="/home/shsoudi/oak/cell-lines/alignment-BWA/dedupped_BAMs_picard/SRC358-N1.sorted.rmdup1.bam"

# Base output directory for Gridss results
base_output_dir="/home/shsoudi/oak/cell-lines/gridss"
# Create the base output directory if it does not exist
mkdir -p "$base_output_dir"

# Sample file containing the sample names or BAM file paths
sample_file="/home/shsoudi/oak/cell-lines/sample.test.txt"


# Loop through each line in the sample file
while IFS= read -r sample_name; do
    # Check if sample_name is a full path to a BAM file or just a sample name
    if [[ "$sample_name" == *".sorted.rmdup1.bam" ]]; then
        bam_file="$sample_name"
        sample_name=$(basename "$bam_file" .sorted.rmdup1.bam)  # Extract sample name from BAM path
    else
        bam_file="$bam_dir/$sample_name.sorted.rmdup1.bam"  # Assume sample_name refers to BAM in bam_dir
    fi
    
    # Check if the BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        echo "BAM file not found for $sample_name: $bam_file"
        continue
    fi
    
    # Create a specific output directory for this sample
    sample_output_dir="$base_output_dir/$sample_name"
    mkdir -p "$sample_output_dir"
    
    # Define the output VCF file path
   output_vcf="$sample_output_dir/${sample_name}_gridss.vcf.gz"
    

    
    # Run Gridss for SV calling
    gridss -r ${reference_genome} -j /home/shsoudi/softwares/gridss/gridss-asset/gridss-2.13.2-gridss-jar-with-dependencies.jar -o "$output_vcf" --skipsoftcliprealignment --threads 8 "$normal_bam" "$bam_file"

        
    
    # Print a message indicating the completion for this sample
    echo "Finished processing $sample_name, results saved in $sample_output_dir"
done < "$sample_file"

# Print a message indicating all samples have been processed
echo "All samples have been processed."