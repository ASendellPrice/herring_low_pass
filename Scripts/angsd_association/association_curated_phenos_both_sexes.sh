#!/bin/bash -l

#SBATCH --array=1-26:1
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -M rackham
#SBATCH -t 5:00:00
#SBATCH -J ANGSD_Association
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools bcftools vcftools/0.1.16

#If directory "angsd_association_both_sex_phenos" does not exist then create it and move into it
#Else 
if [[ ! -d angsd_association_both_sex_phenos ]]
then
    mkdir angsd_association_both_sex_phenos
    cd angsd_association_both_sex_phenos
else
    cd angsd_association_both_sex_phenos
fi

#Determine chromosome name using slurm array task ID
ChrName=chr${SLURM_ARRAY_TASK_ID}

#Create directory for chromosome and move into it
mkdir $ChrName
cd $ChrName

#Set path to VCF file and curated phenotypes file
VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/genotype_likelihoods/${ChrName}/HerringLowPass_GATKMethod_MinMAF0.05_${ChrName}_updatedIDs.vcf.gz
PHENOS=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/curated_phenos_both_sexes.input.txt

#Get basename for VCF and PHENOS
VCF_BASE=$(basename $VCF)
PHENOS_BASE=$(basename $PHENOS)

# Curated phenotypes file for both males and females is structrued as below. Where
# any missing data is denoted by "NA". This file uses the same file format as SNPTEST
# (https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html). This is a 
# space-separated file with two header lines followed by data, as follows:
# Header line: This line lists the name of each column.
# Column type line: This line lists the 'type' of each column:
#     - 0 for the sample ID column
#     - D for a column containing discrete values
#     - P or C for columns containing continuous values (each value must be numerical, or a missing value).
#     - B for a column containing a binary trait (e.g. Male or Female - coded 1 or 2)
# Data lines: Data for each sample. There should be one data line per sample in the genotypes file,
# these must be in the same order that they appear in the genotypes file.

# |-------------------------------------------------------------------------|
# | ID        | sex | age | VS | stomach_weight | body_weight | body_length | <- Header line
# | 0         | B   | D   | D  | P              | P           | P           | <- Column type line
# | 1701-947  | 1   | 6   | 56 | 2.19           | 80.8        | 22.8        | <- Data lines
# | 1701-967  | 1   | 7   | 55 | 1.19           | 68.1        | 20.4        |
# | 1707-1    | 1   | 5   | 57 | 0.70           | 64.5        | 21          |
# | 1707-10   | 1   | 7   | 56 | 0.75           | 70.5        | 21.7        |
# | 1707-100  | 1   | 5   | NA | 0.81           | 67.9        | 21          |
# | 1707-101  | 1   | 5   | 56 | 0.90           | 57          | 18.9        |
# | 1707-102  | 1   | 6   | 56 | 0.46           | 56.9        | 19.2        |
# |-------------------------------------------------------------------------|

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

#Perform association study on body_weight (age and sex = covariates)
ANGSD=/proj/snic2020-2-19/private/herring/users/ashsendell/BIN/angsd/angsd
$ANGSD \
-vcf-gp Subset_${VCF_BASE} \
-fai /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta.fai \
-sampleFile Subset_${PHENOS_BASE} \
-whichPhe body_weight \
-whichCov sex,age \
-doAsso 4 -nInd $N_INDV -doMaf 4 \
-out ${ChrName}_association_body_wieight

#Perform association study on body_length (age and sex = covariates)
$ANGSD \
-vcf-gp Subset_${VCF_BASE} \
-fai /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta.fai \
-sampleFile Subset_${PHENOS_BASE} \
-whichPhe body_length \
-whichCov sex,age \
-doAsso 4 -nInd $N_INDV -minInd 100 -doMaf 4 \
-out ${ChrName}_association_body_length

#Perform association study on VS (age and sex = covariates)
$ANGSD \
-vcf-gp Subset_${VCF_BASE} \
-fai /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta.fai \
-sampleFile Subset_${PHENOS_BASE} \
-whichPhe VS \
-whichCov sex,age \
-doAsso 4 -nInd $N_INDV -doMaf 4 \
-out ${ChrName}_association_VS

#Perform association study on VS (age and sex = covariates)
$ANGSD \
-vcf-gp Subset_${VCF_BASE} \
-fai /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta.fai \
-sampleFile Subset_${PHENOS_BASE} \
-whichPhe stomach_weight \
-whichCov sex,age \
-doAsso 4 -nInd $N_INDV -doMaf 4 \
-out ${ChrName}_association_stomach_weight
