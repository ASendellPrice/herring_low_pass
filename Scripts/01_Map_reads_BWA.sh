#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH --array=1-944:1
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 4:00:00
#SBATCH -J mapping

#Load required modules
module load bioinfo-tools samtools/1.6 bwa

#If directory "mapping" does not exist then create it
if [[ ! -d mapping ]]
then
    mkdir mapping
fi

#Move into mapping directory
cd mapping

#Set paths to reference fasta, sample list and file specifying paths to dirctories containing demultipled reads
REFGENOME=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/sample.list.txt
DEMULTIPLED_READs=/crex/proj/snic2020-2-19/private/herring/users/mats/uppstore2017191/low_pass_de_multiplexing/

#Set sample ID using slurm array task ID
SAMPLE_ID=$(head -n $SLURM_ARRAY_TASK_ID $SAMPLE_LIST | tail -n 1)

R1=${DEMULTIPLED_READs}*/*/*${SAMPLE_ID}_R1.fastq.gz
R2=${DEMULTIPLED_READs}*/*/*${SAMPLE_ID}_R2.fastq.gz
if [ -e $R1 ] && [ -e $R2 ]; then
        bwa mem -t 2 -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID" $REFGENOME $R1 $R2 > ${SAMPLE_ID}.sam  
        samtools view -@ 2 -b -S -o ${SAMPLE_ID}.bam ${SAMPLE_ID}.sam && rm ${SAMPLE_ID}.sam
        samtools sort -@ 2 -o ${SAMPLE_ID}.sort.bam ${SAMPLE_ID}.bam && rm ${SAMPLE_ID}.bam
        samtools index -@ 2 ${SAMPLE_ID}.sort.bam
        samtools stats -@ 2 ${SAMPLE_ID}.sort.bam > ${SAMPLE_ID}.stat
fi
