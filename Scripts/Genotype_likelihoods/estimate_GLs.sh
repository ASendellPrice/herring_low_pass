#!/bin/bash -l

#SBATCH --array=1-26:1
#SBATCH -A snic2022-5-242
#SBATCH -p node -N 1
#SBATCH -M rackham
#SBATCH -t 5-00:00:00
#SBATCH -J Genotype_likelihoods
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION USING ANGSD v.0.930
# Will output the following files per chromosome:
# 1. VCF (with PL field included)
# 2. Beagle (genotype likelihood format)
# 3. MAFs (allele frequencies)

# A. Sendell-Price, Nov 2022
######################################################################################

#STEP 1: Load required modules
ml bioinfo-tools ANGSD/0.921

#STEP 2: Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

#STEP 3: Create directory for chromosome and move into it
mkdir $ChrName
cd $ChrName

#STEP 4: Create bam file list
#Text file containing sample bam paths
ls /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/chrom_bams/${ChrName}/*.sort.bam > chr.bam.list

#STEP 5: Run ANGSD
angsd -b chr.bam.list \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 \
-GL 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doVcf 1 -doGlf 2 -minMaf 0.05 -SNP_pval 1e-6 \
-out HerringLowPass_GATKMethod_MinMAF0.05_${ChrName} -P 20 \
-nThreads 20

# Explanation of above settings:
# ==============================
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -GL = calculate genotype likelihoods (2: using GATK model )
# -doMajorMinor = infer major and minor alleles (1: from GLs)
# -doPost = calculate posterior prob (1: Using frequency as prior)
# -doVcf = output a VCF file (1: yes)
# -doGlf = ouput genotype likelihoods (4: beagle likelihood format)
# -minMaf = minumum minor allele frequency tolerated
# -SNP_pval = significance threshold for determining true polymorphism
# ==============================

#STEP 6: Convert sample names in VCF to something meaningful
#Convert list of bam paths to sample names
for LINE in $(cat chr.bam.list)
do
    basename $LINE | cut -d "." -f 2 >> sample.IDs.txt
done

#Update IDs using bcftools
ml bioinfo-tools bcftools
bcftools reheader --samples sample.IDs.txt \
-o HerringLowPass_GATKMethod_MinMAF0.05_${ChrName}_updatedIDs.vcf.gz \
HerringLowPass_GATKMethod_MinMAF0.05_${ChrName}.vcf.gz

#Remove old VCF file
rm HerringLowPass_GATKMethod_MinMAF0.05_${ChrName}.vcf.gz

