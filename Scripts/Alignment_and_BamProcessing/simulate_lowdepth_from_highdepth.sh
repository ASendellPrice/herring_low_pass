#!/bin/bash -l
#SBATCH -A snic2021-5-8
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M snowy
#SBATCH -t 5-00:00:00
#SBATCH -J DownsampleBams
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools samtools

#Set desired depth of coverage
REQUIRED_DEPTH=1

#For each of the high depth bams do the following ...
for BAM in $(ls /proj/snic2020-2-19/private/herring/alignment/79_individuals/*.bam)
do
#Get ID from bam file name
SAMPLE_ID=$(basename $BAM | cut -d "." -f 1)
#Calculate average depth of coverage for bam
MEAN_DEPTH=$(samtools depth -a $BAM | awk '{sum+=$3} END {print sum/NR}')
#Calculate proportion of reads to keep to acheive desired depth of coverage
PROP_RETAIN=$(awk "BEGIN {print $REQUIRED_DEPTH/$MEAN_DEPTH}")
#Round proportion to 3 decimal places
PROP_RETAIN_3D=$(printf '%.*f\n' 3 $PROP_RETAIN)
#Downsample bam file
samtools view -s $PROP_RETAIN_3D -b $BAM > subsampled_${REQUIRED_DEPTH}X.${SAMPLE_ID}.bam
#Sort bam file
samtools sort subsampled_${REQUIRED_DEPTH}X.${SAMPLE_ID}.bam -o subsampled_${REQUIRED_DEPTH}X.${SAMPLE_ID}.sorted.bam
#Remove un-sorted bam
rm subsampled_${REQUIRED_DEPTH}X.${SAMPLE_ID}.bam
done



