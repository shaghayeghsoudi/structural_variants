#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for add-readgroup-picard
# ---------------------------------------------------------------------
#SBATCH --job-name==picard-addreadgroup
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=40G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=18:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding
#SBATCH --export=ALL
#module purge
# Activate Anaconda work environment for OpenDrift

### NOTE: read groups are necesary for smoove (lumpy)
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""

## Locate the fastqc
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/picard

module load java/17.0.4
module load samtools/1.16.1

## Directories:
# Locate the input data
ROOT_DIR=/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV

# Specify the path to the config file
SCRIPT=/home/users/shsoudi/emoding/envs/picard/share/picard-2.27.5-0
CONFIG=${ROOT_DIR}/sample_ID_config.txt
SAMPLE_ID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIG | cut -f2 -d_)

INPUT_DIR=${ROOT_DIR}/alignment-BWA/dedupped_BAMs_picard


OUTPUT_DIR=/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV/alignment-BWA/dedupped_BAMs_picard/addedRG
mkdir -p ${OUTPUT_DIR}


java -jar ${SCRIPT}/picard.jar AddOrReplaceReadGroups \
      I=${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam \
      O=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1.RG.bam \
      RGID=${SAMPLE_ID} \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=${SAMPLE_ID} \
      RGSM=${SAMPLE_ID} 2> ${SAMPLE_ID}.log
