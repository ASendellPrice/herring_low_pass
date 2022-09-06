#!/bin/bash -l

#SBATCH --array=1-1:1
#SBATCH -A snic2022-5-242
#SBATCH -p node
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 5:00:00
#SBATCH -J whatshapp_phase
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load conda environment
ml conda
source conda_init.sh
#conda create -y -p /home/ashle/conda_envs/impute_genos
conda activate /home/ashle/conda_envs/impute_genos
#conda install whatshap nomkl

#Load other required modules and set path to GLIMPSE directory
ml bioinfo-tools GATK/4.1.1.0 

#Set chromosome name and create / move into chromosome directory
ChrName=chr${SLURM_ARRAY_TASK_ID}

#Create / move into directory "GLIMPSE_imputation"
if [[ ! -d GLIMPSE_imputation ]]
then
    mkdir GLIMPSE_imputation
    cd GLIMPSE_imputation
else
    cd GLIMPSE_imputation
fi

#Make directory for chromosome and a subdirectory called reference
mkdir ${ChrName}
mkdir ${ChrName}/reference

#Move into subdirectory
cd ${ChrName}/reference

#Combine sample GVCFs using GATK CombineGVCFs
REFGENOME=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
GVCFs=$(ls /proj/snic2020-2-19/private/herring/variants/herring_79individuals/*/*.g.vcf.gz | sed -e 's/^/--variant /' | tr "\n" " " )
gatk CombineGVCFs -R $REFGENOME $GVCFs --intervals $ChrName -O ${ChrName}.herring_79individuals.g.vcf.gz

#Call sample genotypes using GATK GenotypeGVCFs
gatk GenotypeGVCFs -R $REFGENOME \
-V ${ChrName}.herring_79individuals.g.vcf.gz -O ${ChrName}.herring_79individuals.vcf.gz

#Remove intermediate GVCF
rm ${ChrName}.herring_79individuals.g.vcf.gz

#NEXT NEED TO APPLY FILTERING TO VCF AND REMOVE INTERMEDIATE FILES
#Apply hard filtering to VCF file
vcftools --gzvcf ${ChrName}.herring_79individuals.vcf.gz \
--remove-indels --min-alleles 2 --max-alleles 2 --mac 2 --minGQ 30 --minDP 10 --max-missing 0.8 \
--recode --out ${ChrName}.herring_79individuals.filtered.vcf

#Remove unfiltered VCF
rm ${ChrName}.herring_79individuals.filtered.vcf

#Phase the filtered VCF using whatshap
BAMS=$(ls /proj/snic2020-2-19/private/herring/alignment/79_individuals/*.bam | tr "\n" " ")
whatshap phase -o ${ChrName}.79_high_depth_samples.phased.vcf \
--reference=${REFGEN} ${ChrName}.79_high_depth_samples.recode.vcf $BAMS

