# How to fix smoove errors

In this post I will show how to solve some Smoove issues I came across. I’m using smoove v0.2.8 installed with conda.

## 1. Duphold error

```
atal.nim(49)            sysFatal
Error: unhandled exception: index -1 not in 0 .. 14130 [IndexDefect]

```

**solution** The duphold version in the conda smoove environment is v0.2.1, need to update it to version 0.2.3, which is not available on conda.
First download duphold from github and check the installation.


```
wget https://github.com/brentp/duphold/releases/download/v0.2.3/duphold
chmod +x ./duphold
./duphold -h # to check if it's installed correctly

```

Add duphold to your path so this version is the version to be used by smoove.

```
PATH=<directory/containing/duphold>:$PATH

```

Since I’m using smoove in a Snakemake pipeline, I add the previous line before my smoove command.

## 2. Svtyper error

```
unrecognized arguments: --max_ci_dist 0

```

**Solution** The svtyper in the smoove conda environment is quite old, you need to update it.

```
conda install svtyper=0.7.0
```

## 3. SR=0 for all variants
When mapping reads with bwa mem, if you use the -M flag, the split reads are maked as secondary and they will not be used by smoove to get split-read support. This will cause all the variants in your VCF file to have SR=0. If you don’t use split-read support in your analysis, you can run smoove as usual, if you want split-read support values in your structural variant file you can use the following solution.
This solution was created by Martijn Derks.

  *It requires python 2, and the pysam and argparse packages as well as samtools*

**solution**

```
module load samtools

sname=`samtools view -H <sample>.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`

python bamgroupreads.py -f -M -i <sample>.bam | samblaster --ignoreUnmated -M -a -e -d <sample>.disc.sam -s <smoove call output dir>/<sample>.split.sam -o /dev/null

grep -v "SAMBLASTER" <smoove call output dir>/<sample>.split.sam > $sname.tmp.sam
mv $sname.tmp.sam <smoove call output dir>/<sample>.split.sam
grep -v "SAMBLASTER" <smoove call output dir>/<sample>.disc.sam > $sname.tmp.sam
mv $sname.tmp.sam <smoove call output dir>/<sample>.disc.sam

samtools sort -@ 12 -O bam <smoove call output dir>/<sample>.split.sam > <smoove call output dir>/$sname.split.bam
samtools sort -@ 12 -O bam <smoove call output dir>/<sample>.disc.sam > <smoove call output dir>/$sname.disc.bam

rm <smoove call output dir>/<sample>.split.sam <smoove call output dir>/<sample>.disc.sam
```

The final split and disc bam files should be in the same directory as the smoove call outdir. Smoove will then use these split and disc bam files for the smoove call step.

## 4.RG error when RG is missing in BAM header

```
load_entry_point('svtyper==0.7.0', 'console_scripts', 'svtyper')()
File "/opt/conda/envs/smoove-env/lib/python2.7/site-packages/svtyper/classic.py", line 575, in cli
sys.exit(main())
File "/opt/conda/envs/smoove-env/lib/python2.7/site-packages/svtyper/classic.py", line 568, in main
args.max_ci_dist)
File "/opt/conda/envs/smoove-env/lib/python2.7/site-packages/svtyper/classic.py", line 150, in sv_genotype
sample = Sample.from_bam(bam_list[i], num_samp, min_lib_prevalence)
File "/opt/conda/envs/smoove-env/lib/python2.7/site-packages/svtyper/parsers.py", line 665, in from_bam
name = bam.header['RG'][0]['SM']
File "pysam/libcalignmentfile.pyx", line 548, in pysam.libcalignmentfile.AlignmentHeader.getitem
KeyError: 'RG'
panic: exit status 1
```

check if your bam/cram files have read-groups. If your BAM file does not have read roup you can use picard tools to add @RG tag

```
java -jar ${SCRIPT}/picard.jar AddOrReplaceReadGroups \
      I=${INPUT_DIR}/${SAMPLE_ID}.rmdup1.bam \
      O=${OUTPUT_DIR}/${SAMPLE_ID}.rmdup1.RG.bam \
      RGID=${SAMPLE_ID} \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=${SAMPLE_ID} \
      RGSM=${SAMPLE_ID} 2> ${SAMPLE_ID}.log
```
