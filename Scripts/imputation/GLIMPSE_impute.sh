#!/bin/bash -l

#SBATCH --array=1-1:1
#SBATCH -A snic2022-5-242
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 1-00:00:00
#SBATCH -J Genotype_Imputation
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

######################################################################################
#STEP 1: Define paths to ANGSD, BAGLE, GLIMPSE and load required modules
######################################################################################

ml bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10 samtools vcftools
BEAGLE4=/proj/snic2020-2-19/private/herring/users/ashsendell/BIN/beagle.27Jan18.7e1.jar
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE

######################################################################################
#STEP 2: Determine chromosome
######################################################################################

ChrName=chr${SLURM_ARRAY_TASK_ID}

######################################################################################
#STEP 3: Create directories
######################################################################################

#If directory "GLIMPSE_imputation" does not exist then create it
if [[ ! -d GLIMPSE_imputation ]]
then
    mkdir GLIMPSE_imputation
fi

#Move into mapping directory
cd GLIMPSE_imputation

#Make directory for focal chrom and move into it
OUT_DIR=$ChrName
mkdir $OUT_DIR
cd $OUT_DIR

######################################################################################
#STEP 4: Create reference panel and list of known SNPs for chromosome
######################################################################################

#Make directory for reference panel and move into it
mkdir reference
cd reference

#Extract polymorphic sites from Mafalda's 91 indv high depth VCF
vcftools --gzvcf /proj/snic2020-2-19/private/herring/users/mafalda/variants/91_indv/vcfs_by_chromosomes_F2/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.vcf.gz \
--min-alleles 2 --max-alleles 2 --remove-indels \
--out herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName} \
--recode

#Phase and index VCF file - this will be used as the reference panel
java -Xmx15000m -jar $BEAGLE4 \
gt=herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.recode.vcf \
out=herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.phased

#Remove un-phased VCF and index phased VCF using tabix
rm herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.recode.vcf
tabix herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.phased.vcf.gz

# Create a list of known variable sites, we will use this later when calculating GLs for low pass data
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' \
herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.phased.vcf.gz \
| bgzip -c > ${ChrName}.sites.tsv.gz
tabix -s1 -b2 -e2 ${ChrName}.sites.tsv.gz -f

######################################################################################
#STEP 5: Split the chromosome into chunks
######################################################################################

#This generates a file containing the imputation chunks and larger chunks including buffers that we will use to run GLIMPSE_phase.
$GLIMPSE_DIR/chunk/bin/GLIMPSE_chunk \
--input herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.phased.vcf.gz \
--region ${ChrName} --window-size 2000000 --buffer-size 200000 \
--output chunks.${ChrName}.txt

#Exit reference directory
cd ../

######################################################################################
#STEP 6: Computing GLs for each individual at variable sites in reference panel
######################################################################################

Make directory for GLs
mkdir genotype_likelihoods
cd genotype_likelihoods

#Create bam file list - a text file containing sample bam paths
for SAMPLE in $(cat ../../resources/imputation/samples_to_impute.txt)
do
    ls ../../chrom_bams/${ChrName}/${ChrName}.${SAMPLE}.sort.bam > ${ChrName}.bam.list
done

#Compute GLs
for BAM in $(cat ${ChrName}.bam.list)
do
samtools index $BAM
SAMPLE_ID=$(basename $BAM | cut -d "." -f 2)
VCF=../reference/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.phased.vcf.gz
TSV=../reference/${ChrName}.sites.tsv.gz
REFGEN=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
OUT=${ChrName}.${SAMPLE_ID}.vcf.gz
bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${ChrName} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}

#Dealing with lack of sample name in bam file
echo $SAMPLE_ID > sample.name
bcftools reheader --samples sample.name -o ${OUT}2 ${OUT}
mv ${OUT}2 ${OUT}

#Index vcf
bcftools index -f ${OUT}
 
#Remove temp file
rm sample.name
done

#Merging genotype likelihoods of multiple individuals
ls ${ChrName}.*.vcf.gz > list.txt
bcftools merge -m none -r ${ChrName} -Oz -o merged.${ChrName}.vcf.gz -l list.txt
tabix merged.${ChrName}.vcf.gz
rm list.txt ${ChrName}.bam.list ${ChrName}.*.gz*

#Exit GLs directory
cd ../

######################################################################################
#STEP 7: Running GLIMPSE
######################################################################################

#To run GLIMPSE_phase we only need to run a job for each imputation chunk. Each job runs on 1 thread.
#Input data are the dataset containing the genotype likelihoods, a reference panel of haplotypes, a genetic map
#and buffered and imputation regions. Here, we run GLIMPSE phase using default parameters that can apply to most datasets. 
#For extremely low coverage and small reference panels, increasing the number of iterations might give more accuracy.

mkdir -p GLIMPSE_impute
cd GLIMPSE_impute

VCF=../genotype_likelihoods/merged.${ChrName}.vcf.gz
REF=../reference/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.variantsOnly.${ChrName}.phased.vcf.gz
MAP=../../../resources/imputation/chrom_maps/${ChrName}.GLIMPSE.gmap.txt

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=${ChrName}.${ID}_GLIMPSE.bcf

#Launch glimpse imputation on chunk in a seperate job, passing the required variables
sbatch ../../scripts/imputation/impute_chunk.sh $VCF $REF $MAP $IRG $ORG $OUT

done < ../reference/chunks.${ChrName}.txt

#Exit impute directory
cd ../

######################################################################################

