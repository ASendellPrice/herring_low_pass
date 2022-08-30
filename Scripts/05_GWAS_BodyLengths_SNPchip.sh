#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p node
#SBATCH -n 1
#SBATCH -M snowy
#SBATCH -t 1-00:00:00
#SBATCH -J GEMMA_BodyLengths
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load conda environment and modules
module load conda
export CONDA_ENVS_PATH=/home/ashle/.conda/envs/
#conda create -n vcf2GWAS
conda activate vcf2GWAS
#conda install vcf2gwas -c conda-forge -c bioconda -c fvogt257
#conda install -c bioconda bcftools 

#Set path to VCF, sample list, PHENOTYPE and COVARIATE files
VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/SNP-chip/VCF/Herring_Array_Hastkar_Filtered_updatedIDs.vcf
SAMPLE_IDS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/SNP-chip/VCF/sample_order_snpchip.txt
PHENOS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/gwas_SNPchip/body_length/bodylength_phenos.csv
COVS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/gwas_SNPchip/body_length/bodylength_covs.csv

#Reorder phenotype and covariate files so their order matches the sample order in the VCF file
head -n 1 $PHENOS > pheno.input
head -n 1 $COVS > covs.input
for SAMPLE in $(cat $SAMPLE_IDS)
do
    grep -w $SAMPLE $PHENOS >> pheno.input
    grep -w $SAMPLE $COVS >> covs.input
done

#Run vcf2gwas
vcf2gwas \
-v $VCF \
-gf /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/SNP-chip/bodylength_gwas_gemma/Clupea_harengus.Ch_v2.0.2.107.gff \
-pf pheno.input \
--topsnp 200 \
-o BodyLength_Age_Sex \
-p 1 -lmm -cf covs.input -c 1 -c 2 






#Set path to VCF, sample list, PHENOTYPE and COVARIATE files
VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/VCFs/AllChromsMerged.phased.imputed.merged.minMAF0.01_NumericIDs.vcf.gz
SAMPLE_IDS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/gwas_GEMMA/sampleIDs_imputed.txt
PHENOS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/gwas_GEMMA/BodyLength_BodyWeight.csv
COVS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/gwas_GEMMA/Sex_Age.csv

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
-gf /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/SNP-chip/bodylength_gwas_gemma/Clupea_harengus.Ch_v2.0.2.107.gff \
-pf pheno.input \
--topsnp 200 \
-o BodyLength_Age_Sex \
-p 1 -lmm -cf covs.input -c 1 -c 2 