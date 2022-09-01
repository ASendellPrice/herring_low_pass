#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH --array=1-26:1
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 1-00:00:00
#SBATCH -J split_bams
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools samtools

#If directory "chrom_bams" does not exist then create it
if [[ ! -d chrom_bams ]]
then
    mkdir chrom_bams
fi

#Move into chrom_bams directory
cd chrom_bams

#Use slurm task ID to define chromosome number
CHROM=chr${SLURM_ARRAY_TASK_ID}

#Make a directory for given chromosome and move into it
mkdir ${CHROM}
cd ${CHROM}

#Set path to file specifying info for 1X data simulated from high pass sequencing
#SampleID   #BamPath
SAMPLE_BAM_INFO=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/simulated_lowpass_info.txt

#Count number of lines in the file
LINE_COUNT=$(cat $SAMPLE_BAM_INFO | wc -l)

#For each sample in the list do the following
for LINE in $(seq 1 $LINE_COUNT)
do
    SAMPLE=$(head -n $LINE $SAMPLE_BAM_INFO | tail -n 1 | cut -f 1)
    BAM_PATH=$(head -n $LINE $SAMPLE_BAM_INFO | tail -n 1 | cut -f 2)
    samtools view -b $BAM_PATH $CHROM > ${CHROM}.${SAMPLE}.sort.bam
    samtools index ${CHROM}.${SAMPLE}.sort.bam
    samtools stats -@ 2 ${CHROM}.${SAMPLE}.sort.bam > ${CHROM}.${SAMPLE}.stat
done
