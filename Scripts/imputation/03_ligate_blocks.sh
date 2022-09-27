#!/bin/bash -l

#SBATCH --array=1-26:1
#SBATCH -A snic2021-5-8
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M snowy
#SBATCH -t 1-00:00:00
#SBATCH -J GLIMPSE_ligate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

######################################################################################

#STEP 1: Define paths to ANGSD, BAGLE, GLIMPSE and load required modules
ml bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE

######################################################################################

#STEP 2: Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

######################################################################################

#STEP 3: Move to chromosome directory
cd $ChrName

######################################################################################

#STEP 4: Ligate multiple chunks together
#The output from the previous step is a VCF/BCF file for each imputed chunk. Now we will merge the chunks into a single file
#for the focal chromosome

mkdir -p GLIMPSE_ligate
cd GLIMPSE_ligate

LST=list.${ChrName}.txt
ls ../GLIMPSE_impute/${ChrName}.*_GLIMPSE.bcf > ${LST}
OUT=${ChrName}.imputed.merged.bcf
${GLIMPSE_DIR}/ligate/bin/GLIMPSE_ligate --input ${LST} --output ${OUT}

######################################################################################

#STEP 5: Phase VCF
IN=${ChrName}.imputed.merged.bcf
OUT=${ChrName}.phased.imputed.merged.vcf
OUT2=${ChrName}.phased.imputed.merged.minMAF0.01.vcf
${GLIMPSE_DIR}/sample/bin/GLIMPSE_sample --input $IN --solve --output $OUT

#STEP 6: Remove non-variant sites and index
bcftools view -q 0.01:minor $OUT | bgzip > ${OUT2}.gz
tabix ${OUT2}.gz

#END


