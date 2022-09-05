#!/bin/bash -l
#SBATCH -A snic2022-5-242
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 2-00:00:00
#SBATCH -J DownsampleBams
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools samtools

#Create directory for simulated low pass data and move into it
if [[ ! -d downsample_high_pass ]]
then
    mkdir downsample_high_pass
fi
cd downsample_high_pass

#Set desired depth of coverage based on average coverage of the 944 low pass samples
TARGET_DEPTH=$(cat ../mapping/*.depth | awk '{sum+=$1} END { print sum/NR}')

#For each of the high depth bams do the following ...
for BAM in $(ls /proj/snic2020-2-19/private/herring/alignment/79_individuals/*.bam)
do
    #Get ID from bam file name
    SAMPLE_ID=$(basename $BAM | cut -d "." -f 1)
    #Calculate average depth of coverage for bam
    MEAN_DEPTH=$(samtools depth -a $BAM | awk '{sum+=$3} END {print sum/NR}')
    #Calculate proportion of reads to keep to acheive desired depth of coverage
    PROP_RETAIN=$(awk "BEGIN {print $TARGET_DEPTH/$MEAN_DEPTH}")
    #Round proportion to 3 decimal places
    PROP_RETAIN_3D=$(printf '%.*f\n' 3 $PROP_RETAIN)
    #Downsample bam file
    samtools view -s $PROP_RETAIN_3D -b $BAM > downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.bam
    #Sort bam file
    samtools sort downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.bam -o downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.sort.bam
    #Remove un-sorted bam
    rm downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.bam
    #Index bam file and output stats
    samtools index -@ 2 downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.sort.bam
    samtools stats -@ 2 downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.sort.bam > downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.stat
    samtools depth downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.sort.bam | awk '{sum+=$3} END { print sum/NR}' > downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.depth
done
