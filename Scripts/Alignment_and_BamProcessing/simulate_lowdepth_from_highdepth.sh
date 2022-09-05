#!/bin/bash -l
#SBATCH -A snic2022-5-242
#SBATCH --array=1-79:1
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 5:00:00
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

#Set high depth bam using slurm array task ID
BAM=$(ls /proj/snic2020-2-19/private/herring/alignment/79_individuals/*.bam | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#Get ID from bam file name
SAMPLE_ID=$(basename $BAM | cut -d "." -f 1)

#Calculate current average depth of bam
REFGENOME=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
REFGENOME_LENGTH=$(cat ${REFGENOME}.fai | awk '{sum+=$2} END {print sum}')
MEAN_DEPTH=$(samtools depth -a $BAM | awk -v len=${REFGENOME_LENGTH} '{sum+=$3} END { print sum/len}')

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

#Calculate sequencing depth of downsampled bam
samtools depth downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.sort.bam \
| awk -v len=${REFGENOME_LENGTH} '{sum+=$3} END { print sum/len}' \
> downsampled_${TARGET_DEPTH}X.${SAMPLE_ID}.depth
