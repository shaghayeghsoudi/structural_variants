#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for Delly
# ---------------------------------------------------------------------
#SBATCH --job-name==Delly
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=70G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=150:00:00
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

export OMP_NUM_THREADS=6
## Locate the fastqc
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/tb-profiler 


module load biology
module load bcftools/1.16

## Directories:
ROOT_DIR=/oak/stanford/groups/emoding/analysis/shaghayegh/shortreads-SV

# Specify the path to the config file
CONFIG=${ROOT_DIR}/sample_ID_config.txt
SAMPLE_ID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIG | cut -f2 -d_)


INPUT_DIR=${ROOT_DIR}/alignment-BWA/dedupped_BAMs_picard
OUT_DIR=${ROOT_DIR}/delly-bcf/${SAMPLE_ID}
mkdir -p ${OUT_DIR}

FASTA=/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa
#delly call -g ${FASTA} ${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam  > ${OUT_DIR}/${SAMPLE_ID}_delly.vcf 



#for name in DEL DUP INV TRA INS; 
#    do delly call -t ${name} -g ${FASTA} ${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam  > ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.vcf
#done

### run delly
for name in DEL DUP INV TRA INS; 
    do delly call -t ${name} -o ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.bcf -g ${FASTA} ${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam  &&
    bcftools view ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.bcf > ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.vcf   #### convert bcf files into vcf

done



## convert bcf too vcf
#for name in DEL DUP INV TRA INS; 
#    bcftools view ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.bcf > ${OUT_DIR}/${SAMPLE_ID}_${name}_delly.vcf
#done 

for file in *vcf;
do bgzip -c ${file} > ${file}.gz && tabix -p vcf ${file}.gz
done

## merge vcf files of SVs for each sample into one vcf
bcftools concat -a *_delly.vcf.gz -o ${SAMPLE_ID}.delly.vcf

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------


