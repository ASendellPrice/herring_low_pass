#!/bin/sh
#SBATCH -A snic2022-5-242
#SBATCH --array=1-1:1
#SBATCH -p node
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 5-00:00:00
#SBATCH -J GenotypeLikelihoods_ANGSD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Run from: /proj/snic2020-2-19/private/herring/users/ashsendell/

#Load required modules
ml bioinfo-tools ANGSD/0.933

#Create directory for analysis and move into it
mkdir GenotypeLikelihoods_ANGSD
cd GenotypeLikelihoods_ANGSD

#Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

#Create directory for chromosome and move into it
OUT_DIR=$ChrName
mkdir $OUT_DIR
cd $OUT_DIR

#Create list of bam files
ls /proj/snic2020-2-19/private/herring/users/ashsendell/Chrom_Bams/${ChrName}/*.sort.bam > bam.list

#Set parameters
REFGENOME=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
BAMLIST=bam.list
OUT=${ChrName}_Ref.is.major_minMAF0.05

#Run ANGSD
angsd -b $BAMLIST -ref $REFGENOME \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -trim 0 \
-doMajorMinor 4 -doMaf 2 -GL 1 -doGlf 2 -SNP_pval 1e-6 -minMaf 0.05 \
-out $OUT -P 5

# Explanation of above settings:
# ==============================
# -ref = reference fasta (used to set Major allele)
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -doMajorMinor = infer major and minor alleles (4: use reference allele as major)
# -doMaf = calculate minor allele frequency (2: fixed major unknown minor)
# -GL = calculate genotype likelihoods (1: using samtools model)
# -doGlf = write GLs to file (2: beagle likelihood file format)
# -SNP_pval = significance threshold for determining true polymorphism
# -minMaf = minumum minor allele frequency tolerated
# ==============================