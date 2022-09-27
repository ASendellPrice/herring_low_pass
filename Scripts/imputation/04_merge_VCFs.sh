#!/bin/bash -l

#SBATCH -A snic2020-5-685
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 10:00:00
#SBATCH -J MergeVCFs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

ml bioinfo-tools bcftools

bcftools concat \
chr*/GLIMPSE_ligate/chr*.phased.imputed.merged.minMAF0.01.vcf.gz \
-o all.chroms.phased.imputed.merged.minMAF0.01.vcf.gz

#Rename chrom so numerical only i.e "1" not "chr1"
zcat all.chroms.phased.imputed.merged.minMAF0.01.vcf.gz | sed 's/^chr//g' | bgzip > temp

mv temp all.chroms.phased.imputed.merged.minMAF0.01.vcf.gz

tabix all.chroms.phased.imputed.merged.minMAF0.01.vcf.gz
