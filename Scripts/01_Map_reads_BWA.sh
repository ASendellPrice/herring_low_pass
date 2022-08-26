#!/bin/bash
#SBATCH -A snic2021-5-8
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 10-00:00:00
#SBATCH -J lane1_BH3G55DRXX
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

module load bioinfo-tools samtools/1.6 bwa


PATH=/crex/proj/snic2020-2-19/private/herring/users/mats/uppstore2017191/low_pass_de_multiplexing/


/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/mapping_lanes.txt


180718_A00181_0044_AH3G5LDRXX   Lane1
180718_A00181_0044_AH3G5LDRXX   Lane2
180718_A00181_0045_BH3G55DRXX   Lane1
180718_A00181_0045_BH3G55DRXX   Lane2


/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/mapping_lanes.txt


list="lane1_BH3G55DRXX.list"
path="/crex/proj/snic2020-2-19/private/herring/users/mats/uppstore2017191/low_pass_de_multiplexing/180718_A00181_0045_BH3G55DRXX/Lane1/"

while read -r line; do
        R1=$path${line}_R1.fastq.gz
        R2=$path${line}_R2.fastq.gz
        if [ -e $R1 ] && [ -e $R2 ]; then
                bwa mem -t 2 -R "@RG\tID:$line\tSM:$line" /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta $R1 $R2 > $line.sam  
                samtools view -@ 2 -b -S -o $line.bam $line.sam && rm $line.sam
                samtools sort -@ 2 -o $line.sort.bam $line.bam && rm $line.bam
                samtools index -@ 2 $line.sort.bam
                samtools stats -@ 2 $line.sort.bam > $line.stat
        fi
done < "$list"
