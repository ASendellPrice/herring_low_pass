#!/bin/bash -l

#SBATCH --array=1-26:1
#SBATCH -A snic2022-5-242
#SBATCH -p core -n 1
#SBATCH -M rackham
#SBATCH -t 1-00:00:00
#SBATCH -J snpEff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

######################################################################################################
# NOTES: HOW TO BUILD A SNPEFF DATABASE
######################################################################################################

#Move into directory containing snpeff
#cd /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/bin/snpEff

#Create directory for new database
#mkdir data/
#mkdir data/Ch_v2.0.2
#cd data/Ch_v2.0.2

#Get annotation files and copy genome
#wget http://ftp.ensembl.org/pub/release-106/gtf/clupea_harengus/Clupea_harengus.Ch_v2.0.2.106.gtf.gz
#zcat Clupea_harengus.Ch_v2.0.2.106.gtf.gz | grep '#!' > genes.gtf
#zcat Clupea_harengus.Ch_v2.0.2.106.gtf.gz | grep -v '#!' | grep -v 'unplaced_' | sed -e 's/^/chr/' >> genes.gtf
#rm Clupea_harengus.Ch_v2.0.2.106.gtf.gz
#ln -s /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta sequences.fa

#Add genome to cofig file
#See: https://pcingola.github.io/SnpEff/se_buildingdb/#add-a-genome-to-the-configuration-file

#Create database:
#cd /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/bin/snpEff
#java -jar snpEff.jar build -gtf22 -v Ch_v2.0.2
######################################################################################################

#Create directory for analysis and move into it
if [[ ! -d snpEff_annotation ]]
then
    mkdir snpEff_annotation
    cd snpEff_annotation
else
    cd snpEff_annotation
fi

#STEP 2: Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

#Set path to snpEff and snpSift jar files and perl file for extracting effects from annotated vcf
snpEff_jar=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/bin/snpEff/snpEff.jar
snpSift_jar=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/bin/snpEff/SnpSift.jar
snpEff_perl=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/bin/snpEff/scripts/vcfEffOnePerLine.pl

#Set path to angsd outputted VCF file
CHROM_VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/genotype_likelihoods/${ChrName}/HerringLowPass_GATKMethod_MinMAF0.01_${ChrName}_updatedIDs.vcf.gz

#Run snpEff (set upstream and downstream interval size to zero using -ud 0)
java -jar $snpEff_jar Ch_v2.0.2 $CHROM_VCF -ud 0 > ${ChrName}.temp.annotated.vcf

#Extract annotations from VCF
#If multiple effects one will be printed per line
cat ${ChrName}.temp.annotated.vcf | $snpEff_perl \
| java -jar $snpSift_jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID:" \
| uniq \
> ${ChrName}_variant_effects.snpEff.txt

#Remove temporary files
rm ${ChrName}.temp.annotated.vcf
