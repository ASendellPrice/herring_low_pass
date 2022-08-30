#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M snowy
#SBATCH -t 1-00:00:00
#SBATCH -J GEMMA_BodyLengths_imputed
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load conda environment and modules
module load conda
export CONDA_ENVS_PATH=/home/ashle/.conda/envs/
#conda create -n vcf2GWAS
conda activate vcf2GWAS
#conda install vcf2gwas -c conda-forge -c bioconda -c fvogt257
#conda install -c bioconda bcftools 


################################################################################################################################################
# RUN ON IMPUTED GENOTYPES
################################################################################################################################################

#Set path to VCF, sample list, PHENOTYPE and COVARIATE files
VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/VCFs/AllChromsMerged.phased.imputed.merged.minMAF0.01_NumericIDs.vcf.gz
SAMPLE_IDS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/gwas_GEMMA/sampleIDs_imputed.txt
PHENOS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/gwas_GEMMA/BodyLength_BodyWeight.csv
COVS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/gwas_GEMMA/Sex_Age.csv
GFF=/proj/snic2020-2-19/private/herring/annotation/Clupea_harengus.Ch_v2.0.2.107.gff

#Reorder phenotype and covariate files so their order matches the sample order in the VCF file
head -n 1 $PHENOS > imputed.bodylength.bodyweight.pheno.input
head -n 1 $COVS > imputed.bodylength.bodyweight.covs.input
for SAMPLE in $(cat $SAMPLE_IDS)
do
    SAMPLE_LINE_NO=$(cat $PHENOS | cut -d "," -f 1 | grep -w $SAMPLE -n | cut -d ":" -f 1)
    if [ ! -z "$SAMPLE_LINE_NO" ]
    then
        head -n $SAMPLE_LINE_NO $PHENOS | tail -n 1 >> imputed.bodylength.bodyweight.pheno.input
        head -n $SAMPLE_LINE_NO $COVS | tail -n 1  >> imputed.bodylength.bodyweight.covs.input
    fi
done

#Run vcf2gwas
vcf2gwas \
-v $VCF \
-gf $GFF \
-pf imputed.bodylength.bodyweight.pheno.input -p 1 -p 2 -p 3 \
--topsnp 200 \
-o BodyLength_Age_Sex_imputed \
-p 1 -lmm -cf imputed.bodylength.bodyweight.covs.input -c 1 -c 2 
