#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH --array=1-26:1
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 2:00:00
#SBATCH -J split_bams
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools samtools

#If directory "chrom_bams" does not exist then create it
if [[ ! -d chrom_bams ]]
then
    mkdir chrom_bams
fi

#Move into mapping directory
cd chrom_bams

#Use slurm task ID to define chromosome number
CHROM=chr${SLURM_ARRAY_TASK_ID}

#Make a directory for given chromosome and move into it
mkdir ${CHROM}
cd ${CHROM}

#Set path to sample list and directory containing sample BAMs
SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/sample.list.txt
BAM_DIRECTORY=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/mapping

while read -r line; do
    samtools view -b ${BAM_DIRECTORY}/${line}.sort.bam ${CHROM} > chr${CHROM}.${line}.sort.bam
    samtools index chr${CHROM}.${line}.sort.bam
    samtools stats -@ 2 chr${CHROM}.${line}.sort.bam > chr${CHROM}.${SAMPLE_ID}.stat
done < "$SAMPLE_LIST"
