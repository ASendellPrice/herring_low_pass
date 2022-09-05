#!/usr/bin/env Rscript
setwd("/Users/ashse759/Dropbox/Uppsala_PostDoc/herring/")

library(ggplot2)
library(ggthemes)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(plotly)
library(factoextra)
library(xlsx)
library(ggthemr)
ggthemr('flat', layout = 'clean')


#Load bodylength GWAS results for imputed genotypes
imputed_pwald <- read.delim("BodyLength_Age_Sex_VS_imputed/output/lmm/body_length/body_length_20220830_125617/body_length_mod_sub_imputed.bodylength.bodyweight.vs.pheno.input.part1_AllChromsMerged.phased.imputed.merged.minMAF0.01_NumericIDs.assoc.txt") %>%
  select(chr, ps, p_wald) %>% rename(p_wald_imputed = p_wald)

#Load bodylength GWAS results for imputed genotypes
SNPchip_pwald <- read.delim("BodyLength_Age_Sex_VS_snpchip/output/lmm/body_length/body_length_20220830_115130/body_length_mod_sub_snpchip.bodylength.bodyweight.vs.pheno.input.part1_Herring_Array_Hastkar_Filtered_NumericIDs.assoc.txt") %>%
  select(chr, ps, p_wald) %>% rename(p_wald_snpchip = p_wald)

#Combine the two datasets
combined <- full_join(imputed_pwald, SNPchip_pwald, by = c("chr", "ps"))
rm(imputed_pwald, SNPchip_pwald)

#Create new column specifying postion of snp in dataframe from 1 to n
combined$n <- 1:nrow(combined)

#Split up the dataframes again
imputed_pwald <- na.omit(combined %>% select(chr, ps, n, p_wald_imputed) %>% rename(p_wald = p_wald_imputed))
SNPchip_pwald <- na.omit(combined %>% select(chr, ps, n, p_wald_snpchip) %>% rename(p_wald = p_wald_snpchip))

#Add column specifying the dataset
imputed_pwald$dataset <- "imputed"
SNPchip_pwald$dataset <- "snp_chip"

#Combine once again
combined <- rbind(imputed_pwald, SNPchip_pwald)
rm(imputed_pwald, SNPchip_pwald)

#Calculate -log10 of pval
combined$negLog10 <- log10(combined$p_wald) * -1

#Load list of significantly associated SNPs
imputed_top_associated <- read.csv("BodyLength_Age_Sex_VS_imputed/output/lmm/body_length/body_length_20220830_125617/best_p-values/p_wald_body_length_mod_sub_imputed.bodylength.bodyweight.vs.pheno.input.part1_AllChromsMerged.phased.imputed.merged.minMAF0.01_NumericIDs_top001%.csv") %>%
  select(chr, ps)

#Add info from combined dataframe to significantly associated SNP list
imputed_top_associated <- left_join(imputed_top_associated, combined)

#For each chromosome find centre in n - used to add x axis labels
axisdf <- combined %>% filter(dataset == "imputed") %>% group_by(chr) %>% summarize(center=(max(n) + min(n) ) / 2 )

#Set up whole genome plot (imputed)
P_imputed <- ggplot(filter(combined, dataset == "imputed"), aes(x=n, y=negLog10, color=as.factor(chr))) +
  
  # Show all points
  geom_point(size = 0.3) +
  scale_color_manual(values = rep(c("#242565","#46557E"), 26 )) +
  #geom_hline(yintercept = min(GEMMA.sig.SNPs.info$negLog10), linetype = "dashed", colour = "grey37") +
  geom_point(data = filter(imputed_top_associated, dataset == "imputed"), colour = "#F71927", size = 0.3) +
  #geom_text(data = clusterMid, inherit.aes = FALSE, aes(x = center, y = ypos, label = paste0("#",clusterID)), size = 3, colour = "black") +
  #geom_label_repel(data = na.omit(GEMMA.sig.SNPs.info), aes(label = gene, fill = effect), force = 15, nudge_y = 0, size = 2.5, na.rm = TRUE, max.overlaps = 50, segment.colour = "black", colour = "white", fontface = 'bold', alpha = 0.75) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center , expand = c(0.01, 0.01)) +
  labs(x = "Chromosome", y = "-Log10(p-value)") +
  
  # Customise the theme:
  theme(legend.position="none", panel.border = element_blank())

#Save plot as png file
ggsave(P_imputed, file="BodyLength_Age_Sex_VS_imputed/output/lmm/body_length/body_length_20220830_125617/gwas_bodylength_imputed.png", width = 25, height = 6, units = "cm")


#Set up whole genome plot (snpchip)
snp_chip <- combined %>% filter(dataset == "snp_chip") %>% filter(chr > 0)
P_snpchip <- ggplot(snp_chip, aes(x=n, y=negLog10, color=as.factor(chr))) +
  
  # Show all points
  geom_point(size = 0.3) +
  scale_color_manual(values = rep(c("#242565","#46557E"), 26 )) +
  #geom_hline(yintercept = min(GEMMA.sig.SNPs.info$negLog10), linetype = "dashed", colour = "grey37") +
  geom_point(data = filter(imputed_top_associated, dataset == "snp_chip"), colour = "#F71927", size = 0.3) +
  #geom_text(data = clusterMid, inherit.aes = FALSE, aes(x = center, y = ypos, label = paste0("#",clusterID)), size = 3, colour = "black") +
  #geom_label_repel(data = na.omit(GEMMA.sig.SNPs.info), aes(label = gene, fill = effect), force = 15, nudge_y = 0, size = 2.5, na.rm = TRUE, max.overlaps = 50, segment.colour = "black", colour = "white", fontface = 'bold', alpha = 0.75) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center , expand = c(0.01, 0.01)) +
  labs(x = "Chromosome", y = "-Log10(p-value)") +
  
  # Customise the theme:
  theme(legend.position="none", panel.border = element_blank())

#Save plot as png file
ggsave(P_snpchip, file="BodyLength_Age_Sex_VS_snpchip/output/lmm/body_length/body_length_20220830_115130/gwas_bodylength_snpchip.png", width = 25, height = 6, units = "cm")

