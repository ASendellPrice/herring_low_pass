#!/bin/bash -l

#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -M rackham
#SBATCH -t 5:00:00
#SBATCH -J PCA_tshr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools ANGSD/0.921 PCAngsd/0.982

#
mkdir PCA_tshr_region
cd PCA_tshr_region

#Create bam file list
#Text file containing sample bam paths
ls /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/chrom_bams/chr15/*.sort.bam > chr.bam.list

#STEP 5: Run ANGSD
angsd -b chr.bam.list \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 \
-GL 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doVcf 1 -doGlf 2 -minMaf 0.05 -SNP_pval 1e-6 \
-r chr15:8850000-8950000 -out HerringLowPass_GATKMethod_MinMAF0.05_TSHR_region -P 20 \
-nThreads 20

#Convert list of bam paths to sample names
for LINE in $(cat chr.bam.list)
do
    basename $LINE | cut -d "." -f 2 >> sample.IDs.txt
done

#Run PCAngsd
pcangsd.py -beagle HerringLowPass_GATKMethod_MinMAF0.05_TSHR_region.beagle.gz \
-iter 10000 -o HerringLowPass_GATKMethod_MinMAF0.05_TSHR_region

