#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for Lumpy
# ---------------------------------------------------------------------
#SBATCH --job-name==Lumpy-Smoove
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=60G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=148:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding
#SBATCH --export=ALL
#module purge
# Activate Anaconda work environment for OpenDrift


# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""

#module load system openssl 
## Locate the smoove conda
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/smoove-only

## Directories:
ROOT_DIR=/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV
BED=/oak/stanford/groups/emoding/analysis/shaghayegh/resources/lumpy/ceph18.b37.lumpy.exclude.2014-01-15.bed
# Specify the path to the config file
CONFIG=${ROOT_DIR}/sample_ID_config.txt
SAMPLE_ID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIG | cut -f2 -d_)


#INPUT_DIR=${ROOT_DIR}/alignment-BWA
INPUT_DIR=${ROOT_DIR}/alignment-BWA/dedupped_BAMs_picard/addedRG

OUT_DIR=${ROOT_DIR}/results-smoove/${SAMPLE_ID}
mkdir -p ${OUT_DIR}

LOG_DIR=${ROOT_DIR}/results-smoove/${SAMPLE_ID}/stderr
mkdir -p ${LOG_DIR}

FASTA=/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa


smoove call -o ${OUT_DIR} --exclude ${BED} --name ${SAMPLE_ID} --fasta ${FASTA} -p 10 ${INPUT_DIR}/${SAMPLE_ID}.rmdup1.RG.bam 2> ${LOG_DIR}/${SAMPLE_ID}.log

#smoove call --outdir test_mansi/SRR6745447_output/smoove_output/ --name SRR6745447 --fasta part_II/references/hg38.fa -p 1 --genotype part_II/output/fq2bam_no_indels/SRR6745447_sorted.bam

##--genotype is causing issue because of svtyper

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------


