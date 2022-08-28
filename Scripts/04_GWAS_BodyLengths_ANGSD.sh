#!/bin/sh
#SBATCH -A snic2022-5-242
#SBATCH --array=1-1:1
#SBATCH -p node
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 1-00:00:00
#SBATCH -J gwas_angsd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools ANGSD/0.933

#If directory "gwas_bodylength_angsd" does not exist then create it
if [[ ! -d gwas_angsd ]]
then
    mkdir gwas_bodylength_angsd
fi

#Move into gwas_bodylength_angsd directory
cd gwas_bodylength_angsd

#Determine chromosome name using slurm array task id
CHROM=chr${SLURM_ARRAY_TASK_ID}

#Create directory for chromosome and move into it
mkdir $CHROM
cd $CHROM

#Set path to chromosome bam directory
BAM_DIRECTORY=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/chrom_bams/${CHROM}

#Set path to reference fasta
REFGENOME=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

#Set path to files specifying samples to use, their phenotypic measurements and covariates
SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/body_lengths_samples.txt
PHENOTYPE=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/body_lengths.txt
COVARIATES=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/body_lengths_sex_age.txt

#Set output prefix
OUT=Association_BodyLength_${CHROM}

#Create a bam input list by looping through the sample list
#Create list of chromosome bam files
while read -r line; do
    ls ${BAM_DIRECTORY}/${CHROM}.*-${line}.sort.bam >> bam.list
done < "$SAMPLE_LIST"

#Run ANGSD
angsd -b bam.list -ref $REFGENOME \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -trim 0 \
-doMajorMinor 4 -doMaf 2 -GL 1 -doGlf 2 -SNP_pval 1e-6 -minMaf 0.05 \
-P 20 -nThreads 20 -out $OUT 


#\
#-doAsso 2 -yQuant $PHENOTYPE -cov $COVARIATES \
#-out $OUT -P 5
