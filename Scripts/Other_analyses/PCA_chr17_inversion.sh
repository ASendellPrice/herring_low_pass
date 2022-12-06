#!/bin/bash -l

#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -M rackham
#SBATCH -t 5:00:00
#SBATCH -J PCA_chr17
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools ANGSD/0.921 PCAngsd/0.982

#Create directory for analysis
mkdir PCA_chr17_inversion
cd PCA_chr17_inversion

#Create bam file list
#Text file containing sample bam paths
for SAMPLE in $(cat ../resources/samples_for_spawing_PCA.txt)
do
ls /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/chrom_bams/chr17/*.${SAMPLE}.sort.bam >> bam.list
done

#Create beagle file for region of interest
angsd -b  bam.list \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 \
-GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 \
-r chr17:25800000-27500000 -out HerringLowPass_GATKMethod_chr17_inversion

#Create covariance matrix
pcangsd.py -beagle HerringLowPass_GATKMethod_chr17_inversion.beagle.gz \
-minMaf 0.01 -iter 10000 -o HerringLowPass_GATKMethod_chr17_inversion

#Convert bam.list to list of sample names
for LINE in $(cat bam.list)
do
    basename $LINE | cut -d "." -f 2 >> sample.IDs.txt
done
