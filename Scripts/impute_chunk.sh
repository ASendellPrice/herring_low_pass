#!/bin/bash -l

#SBATCH -A snic2021-5-8
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M snowy
#SBATCH -t 10-00:00:00
#SBATCH -J Genotype_Imputation_Chunk
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

######################################################################################
#STEP 1: Define paths to ANGSD, BAGLE, GLIMPSE and load required modules
######################################################################################

ml bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10 samtools vcftools
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE

######################################################################################
#STEP 2: Retrieve variables from batch submission
######################################################################################

VCF=${1}
REF=${2}
MAP=${3}
IRG=${4}
ORG=${5}
OUT=${6}

######################################################################################
#STEP 3: Run GLIMPSE on chunk
######################################################################################

${GLIMPSE_DIR}/phase/bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT} \
--burnin 50 --main 15 
bcftools index -f ${OUT}
