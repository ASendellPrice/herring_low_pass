#!/bin/sh
#SBATCH -A snic2022-5-242
#SBATCH --array=1-1:1
#SBATCH -p node
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 5-00:00:00
#SBATCH -J GenotypeLikelihoods_ANGSD
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools ANGSD/0.933

#If directory "genotype_likelihoods" does not exist then create it
if [[ ! -d genotype_likelihoods ]]
then
    mkdir genotype_likelihoods
fi

#Move into chrom_bams directory
cd genotype_likelihoods

#Determine chromosome name using slurm array task id
$CHROM=chr${SLURM_ARRAY_TASK_ID}

#Create directory for chromosome and move into it
mkdir $CHROM
cd $CHROM

#Set path to sample list and chromosome bam directory
SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/sample.list.txt
BAM_DIRECTORY=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/chrom_bams/${CHROM}

#Create list of chromosome bam files
while read -r line; do
    ls ${BAM_DIRECTORY}/${CHROM}.${line}.sort.bam >> bam.list
done < "$SAMPLE_LIST"

#Set parameters
REFGENOME=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
OUT=${CHROM}_Ref.is.major_minMAF0.05

#Run ANGSD
angsd -b bam.list -ref $REFGENOME \
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