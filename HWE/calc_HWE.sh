

#Load required modules
ml bioinfo-tools vcftools/0.1.16

#Set path to SNPchip VCF file
VCF=/proj/snic2020-2-19/private/herring/users/ashsendell/herring_low_pass/VCFs/Herring_Array_Hastkar_Filtered_NumericIDs.vcf

#Extract TSHR region from VCF file
vcftools --vcf $VCF --chr 15 --from-bp 8500000 --to-bp 9500000 --recode --stdout \
> Herring_Array_Hastkar_Filtered_NumericIDs_TSHRregion.vcf

#Chr6 inversion
vcftools --vcf $VCF --chr 6 --from-bp 22200000 --to-bp 24800000 --recode --stdout \
> Herring_Array_Hastkar_Filtered_NumericIDs_chr6inversion.vcf

#Chr12 inversion
vcftools --vcf $VCF --chr 12 --from-bp 17800000 --to-bp 25600000 --recode --stdout \
> Herring_Array_Hastkar_Filtered_NumericIDs_chr12inversion.vcf

#Chr17 inversion
vcftools --vcf $VCF --chr 17 --from-bp 25800000 --to-bp 27500000 --recode --stdout \
> Herring_Array_Hastkar_Filtered_NumericIDs_chr17inversion.vcf

#Chr23 inversion
vcftools --vcf $VCF --chr 23 --from-bp 16300000 --to-bp 17500000 --recode --stdout \
> Herring_Array_Hastkar_Filtered_NumericIDs_chr23inversion.vcf

