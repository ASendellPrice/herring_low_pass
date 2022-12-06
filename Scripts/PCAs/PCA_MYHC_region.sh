#!/bin/bash -l

#SBATCH -A snic2022-5-241
#SBATCH -p core -n 1
#SBATCH -M rackham
#SBATCH -t 14:00:00
#SBATCH -J PCA_MYHC_region
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
ml bioinfo-tools ANGSD/0.921 PCAngsd/0.982

#Create directory for analysis
if [[ ! -d pcangsd_output ]]
then
    mkdir pcangsd_output
    cd pcangsd_output
else
    cd pcangsd_output
fi

#Set path to file specifying samples to include and reference genome
SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/resources/samples_for_spawing_PCA.txt
REF=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

#Set chromosome name, output prefix and path to site list
CHROM=chr12
OUT=HerringLowPass_GATKMethod_MYHC_region
START=15763753
END=16047860

#Create bam file list
#Text file containing sample bam paths
for SAMPLE in $(cat $SAMPLE_LIST)
do
ls /proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/chrom_bams/${CHROM}/*.${SAMPLE}.sort.bam >> ${OUT}.bam.list
done

#Create beagle file for region of interest
angsd -b ${OUT}.bam.list -ref $REF \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 \
-GL 2 -doMajorMinor 4 -doMaf 1 -minMaf 0.01 -doGlf 2 \
-r ${CHROM}:${START}-${END} -out $OUT

#Create covariance matrix
pcangsd.py -beagle ${OUT}.beagle.gz \
-minMaf 0.01 -iter 10000 -o $OUT

#Convert bam list to list of sample names
for LINE in $(cat ${OUT}.bam.list)
do
    basename $LINE | cut -d "." -f 2 | sed 's/1701-/Hastskar_/g' | sed 's/1707-/Hastskar_/g' >> ${OUT}.sample.IDs.txt
done
