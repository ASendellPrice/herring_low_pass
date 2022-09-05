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

#For each of the low pass bam files do the the following:
for BAM in $(ls ../../mapping/*.sort.bam)
do
    BASE_NAME=$(basename $BAM)
    samtools view -b $BAM ${CHROM} > ${CHROM}.${BASE_NAME}
    samtools index ${CHROM}.${BASE_NAME}
done

#Now do the same for the simulated low pass bams ...
for BAM in $(ls ../../simulated_low_pass/*.sort.bam)
do
    BASE_NAME=$(basename $BAM)
    samtools view -b $BAM ${CHROM} > ${CHROM}.${BASE_NAME}
    samtools index ${CHROM}.${BASE_NAME}
done
