#!/bin/bash -l

#SBATCH --array=19-19:1
#SBATCH -A snic2022-5-242
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 10:00:00
#SBATCH -J Genotype_Imputation
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

######################################################################################
#Define paths to GLIMPSE and load required modules
######################################################################################

#Load other required modules and set path to GLIMPSE directory
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE
ml bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10 samtools vcftools

######################################################################################
#Determine chromosome using slurm array job id and move into directory
######################################################################################

ChrName=chr${SLURM_ARRAY_TASK_ID}
cd GLIMPSE_imputation/${ChrName}
DIRECTORY=$(pwd)

######################################################################################
# Split the chromosome into chunks
######################################################################################

#Move into referencd directory
cd reference

#Index VCF (this is the reference panel)
bcftools index ${ChrName}.herring_79individuals.filtered.phased.vcf.gz

# Create a list of known variable sites, we will use this later when calculating GLs for low pass data
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' \
${ChrName}.herring_79individuals.filtered.phased.vcf.gz \
| bgzip -c > ${ChrName}.sites.tsv.gz
tabix -s1 -b2 -e2 ${ChrName}.sites.tsv.gz -f

#Generate a file containing the imputation chunks and larger chunks including buffers that we will use to run GLIMPSE_phase.
$GLIMPSE_DIR/chunk/bin/GLIMPSE_chunk \
--input ${ChrName}.herring_79individuals.filtered.phased.vcf.gz \
--region ${ChrName} --window-size 2000000 --buffer-size 200000 \
--output chunks.${ChrName}.txt

#Exit reference directory
cd ../

######################################################################################
#Computing GLs for each individual at variable sites in reference panel
######################################################################################

#Make directory for GLs
mkdir genotype_likelihoods
cd genotype_likelihoods

#Set path to reference assembly
REFGEN=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

#Compute GLs from low pass bams
for BAM in $(ls ../../../chrom_bams/${ChrName}/*.bam)
do
    SAMPLE_ID=$(basename $BAM | sed "s/${ChrName}.//g" | sed "s/.sort.bam//g")
    VCF=../reference/${ChrName}.herring_79individuals.filtered.phased.vcf.gz
    TSV=../reference/${ChrName}.sites.tsv.gz
    OUT=${ChrName}.${SAMPLE_ID}.vcf.gz
    bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${ChrName} -q 20 -Q 20 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}

    #Dealing with lack of sample name in bam file
    echo $SAMPLE_ID > sample.name
    bcftools reheader --samples sample.name -o ${OUT}2 ${OUT}
    rm sample.name
    mv ${OUT}2 ${OUT}

    #Index vcf
    bcftools index -f ${OUT}
done

#Merging genotype likelihoods of multiple individuals
ls ${ChrName}.*.vcf.gz > list.txt
bcftools merge -m none -r ${ChrName} -Oz -o merged.${ChrName}.vcf.gz -l list.txt
tabix merged.${ChrName}.vcf.gz
rm list.txt ${ChrName}.*.gz*

#Exit GLs directory
cd ../

######################################################################################
#Running GLIMPSE
######################################################################################

#To run GLIMPSE_phase we only need to run a job for each imputation chunk. Each job runs on 1 thread.
#Input data are the dataset containing the genotype likelihoods, a reference panel of haplotypes, a genetic map
#and buffered and imputation regions. Here, we run GLIMPSE phase using default parameters that can apply to most datasets. 
#For extremely low coverage and small reference panels, increasing the number of iterations might give more accuracy.

mkdir -p GLIMPSE_impute
cd GLIMPSE_impute

VCF=${DIRECTORY}/genotype_likelihoods/merged.${ChrName}.vcf.gz
REF=${DIRECTORY}/reference/${ChrName}.herring_79individuals.filtered.phased.vcf.gz
MAP=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/imputation/chrom_maps/${ChrName}.GLIMPSE.gmap.txt

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=${ChrName}.${ID}_GLIMPSE.bcf

#Launch glimpse imputation on chunk in a seperate job, passing the required variables
sbatch ../../../Scripts/imputation/impute_chunk.sh $VCF $REF $MAP $IRG $ORG $OUT

done < ../reference/chunks.${ChrName}.txt

######################################################################################




for CHROM in $(seq 2 3)
do

ChrName=chr${CHROM}
cd $ChrName
DIRECTORY=$(pwd)
cd GLIMPSE_impute

VCF=${DIRECTORY}/genotype_likelihoods/merged.${ChrName}.vcf.gz
REF=${DIRECTORY}/reference/${ChrName}.herring_79individuals.filtered.phased.vcf.gz
MAP=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/imputation/chrom_maps/${ChrName}.GLIMPSE.gmap.txt

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=${ChrName}.${ID}_GLIMPSE.bcf

#Launch glimpse imputation on chunk in a seperate job, passing the required variables
sbatch ../../../Scripts/imputation/impute_chunk.sh $VCF $REF $MAP $IRG $ORG $OUT

done < ../reference/chunks.${ChrName}.txt

cd ../../

done