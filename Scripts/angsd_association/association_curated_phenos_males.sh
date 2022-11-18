#!/bin/bash -l

#SBATCH --array=1-26:1
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -M rackham
#SBATCH -t 5-00:00:00
#SBATCH -J ANGSD_Association
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools bcftools vcftools/0.1.16

#If directory "angsd_association_both_sex_phenos" does not exist then create it and move into it
#Else 
if [[ ! -d angsd_association_male_phenos ]]
then
    mkdir angsd_association_male_phenos
    cd angsd_association_male_phenos
else
    cd angsd_association_male_phenos
fi

#Determine chromosome name using slurm array task ID
ChrName=chr${SLURM_ARRAY_TASK_ID}

#Create directory for chromosome and move into it
mkdir $ChrName
cd $ChrName

#Set path to VCF file and curated phenotypes file
VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/genotype_likelihoods/${ChrName}/HerringLowPass_GATKMethod_MinMAF0.05_${ChrName}_updatedIDs.vcf.gz
PHENOS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/curated_phenos_males.input.txt

#Get basename for VCF and PHENOS
VCF_BASE=$(basename $VCF)
PHENOS_BASE=$(basename $PHENOS)

#From curated phenotypes file extract list of sample IDs:
tail -n +3 $PHENOS | cut -f 1 > phenotyped_samples.IDs.txt

#Filter VCF file keeping only phenotyped samples
zcat $VCF \
| sed 's/VCFv4.2(angsd version)/VCFv4.2/g' \
| vcftools --vcf - --keep phenotyped_samples.IDs.txt --recode --stdout \
| bgzip > Subset_${VCF_BASE} # <- can angsd pass gzipped vcfs???
tabix Subset_${VCF_BASE}

#Get sample order in the VCF file
bcftools query -l Subset_${VCF_BASE} \
> sample_order_in_subsetted_vcf.txt

#Count number of samples in that list (needed later)
N_INDV=$(cat sample_order_in_subsetted_vcf.txt | wc -l)

#Reorder curated phenotypes file so that samples are in the same order that they
#appear in the VCF file.
head -n 2 $PHENOS > Subset_${PHENOS_BASE}
for SAMPLE in $(cat sample_order_in_subsetted_vcf.txt)
do
    cat $PHENOS | grep -w $SAMPLE >> Subset_${PHENOS_BASE}
done

#Perform association study on gonad_weight (age = covariate)
ANGSD=/proj/snic2020-2-19/private/herring/users/ashsendell/BIN/angsd/angsd
$ANGSD \
-vcf-gp Subset_${VCF_BASE} \
-fai /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta.fai \
-sampleFile Subset_${PHENOS_BASE} \
-whichPhe gonad_weight \
-whichCov age \
-doAsso 4 -nInd $N_INDV -doMaf 4 \
-out ${ChrName}_association_male_gonad_weight

