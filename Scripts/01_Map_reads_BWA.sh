#!/bin/bash
#SBATCH -A snic2021-5-8
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 10-00:00:00
#SBATCH -J mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
module load bioinfo-tools samtools/1.6 bwa

#Set path to file specifying directories containing demultiplexed reads
PATHS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/read_directories.txt

#Set directory path according to slurm array ID
READ_DIRECTORY=$(head -n $SLURM_ARRAY_TASK_ID $PATHS | tail -n 1)

#Extract read ID (e.g. 180718_A00181_0044_AH3G5LDRXX) and lane number from path
READ_ID=$(echo $READ_DIRECTORY | cut -d "/" -f 11)
LANE=$(echo $READ_DIRECTORY | cut -d "/" -f 12)

#Based on the files in the read directory create a list of file names
ls $READ_DIRECTORY | grep "R1" | cut -d "_" -f 1 > ${READ_ID}.${LANE}.sample.list

#Set path to reference fasta
REFGENOME=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

while read -r line; do
        R1=${READ_DIRECTORY}${line}_R1.fastq.gz
        R2=${READ_DIRECTORY}${line}_R1.fastq.gz
        if [ -e $R1 ] && [ -e $R2 ]; then
                bwa mem -t 2 -R "@RG\tID:$line\tSM:$line" $REFGENOME $R1 $R2 > ${line}_${READ_ID}_${LANE}.sam  
                samtools view -@ 2 -b -S -o ${line}_${READ_ID}_${LANE}.bam ${line}_${READ_ID}_${LANE}.sam && rm ${line}_${READ_ID}_${LANE}.sam
                samtools sort -@ 2 -o ${line}_${READ_ID}_${LANE}.sort.bam ${line}_${READ_ID}_${LANE}.bam && rm ${line}_${READ_ID}_${LANE}.bam
                samtools index -@ 2 ${line}_${READ_ID}_${LANE}.sort.bam
                samtools stats -@ 2 ${line}_${READ_ID}_${LANE}.sort.bam > ${line}_${READ_ID}_${LANE}.stat
        fi
done < "${READ_ID}.${LANE}.sample.list"
