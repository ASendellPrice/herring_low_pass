#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH --array=1-4:1
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 5-00:00:00
#SBATCH -J mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashley.sendell-price@imbim.uu.se

#Load required modules
module load bioinfo-tools samtools/1.6

#Move into mapping directory
cd mapping


samtools stats -@ 2