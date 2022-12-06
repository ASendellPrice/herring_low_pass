setwd("/Users/ashse759/Dropbox/Uppsala_PostDoc/herring_low_pass/PCA_tshr_region/")

library(dplyr)
library(ggplot2)
library(ggthemes)
library(e1071)

#Read in covariance matrix for tshr region and compute eigenvalues and eigenvectors
cov <- as.matrix(read.table("HerringLowPass_GATKMethod_TSHR_region.cov"))
e <- eigen(cov)

#Extract first two PCs
PCs <- as.data.frame(e$vectors[,1:2])

#Add sample IDs to dataframe
PCs$sample_ID <- readLines("sample_IDs.txt") 

#Add spawning season info (extract from sample name if present)
PCs$spawning.season <- NA

for (SAMPLE in 1:nrow(PCs)){
  if (stringr::str_detect(PCs$sample_ID[SAMPLE], 'Spring') == TRUE){
    PCs$spawning.season[SAMPLE] <- "Spring" 
  }
  else if (stringr::str_detect(PCs$sample_ID[SAMPLE], 'Autumn') == TRUE){
    PCs$spawning.season[SAMPLE] <- "Autumn" 
  } 
  else if (stringr::str_detect(PCs$sample_ID[SAMPLE], 'Summer') == TRUE){
    PCs$spawning.season[SAMPLE] <- "Summer" 
  } 
  else if (stringr::str_detect(PCs$sample_ID[SAMPLE], 'Winter') == TRUE){
    PCs$spawning.season[SAMPLE] <- "Winter" 
  } 
}

#Plot PC1 and 2
#Based on this Autumn and spring are segregating along PC1,
#though there are many of the hästskär samples are intermediate (the NAs)
P1 <- ggplot(data = PCs, aes(x = V1, y = V2, colour = spawning.season)) +
  geom_point() +
  labs(x = paste0("PC1 (",round(e$values[1],2),"%)"), 
       y = paste0("PC2 (",round(e$values[2],2),"%)"), 
       title = "TSHR region: 8.85-8.95Mb") +
  scale_color_manual(values = c("#3498db","#2ecc71","#f1c40f")) +
  theme_bw()
P1

#Assign samples to 1 of 3 clusters along PC1 axis
centers <- cmeans(PCs$V1, 3)$centers
centers <- sort(centers)
cm <- cmeans(PCs$V1, centers = centers)

#Create data.frame storing cluster membership probabilities plus most likely assignment
membership <- as.data.frame(cm$membership)
membership$BestCluster <- cm$cluster

#Assign samples to a cluster so long as its assignment probability is at least 95%
PCs$cluster <- NA
for (SAMPLE in 1:nrow(membership)){
  cluster_assigned <- membership$BestCluster[SAMPLE]
  if (membership[SAMPLE,cluster_assigned] >= 0.80){
    PCs$cluster[SAMPLE] <- membership$BestCluster[SAMPLE]
  }
}

#Update cluster names
PCs$spawning.season.pc.based <- NA
for (SAMPLE in 1:nrow(PCs)){
  if (!is.na(PCs$cluster[SAMPLE]) && PCs$cluster[SAMPLE] == 1){
    PCs$spawning.season.pc.based[SAMPLE] <- "Spring-like"
  } 
  else if (!is.na(PCs$cluster[SAMPLE]) && PCs$cluster[SAMPLE] == 2){
    PCs$spawning.season.pc.based[SAMPLE] <- "Intermediate"
  }
  else if (!is.na(PCs$cluster[SAMPLE]) && PCs$cluster[SAMPLE] == 3){
    PCs$spawning.season.pc.based[SAMPLE] <- "Autumn-like"
  }
}
  
#Plot PCs this time colouring samples by infered spawing type
P2 <- ggplot(data = PCs, aes(x = V1, y = V2, colour = spawning.season.pc.based)) +
  geom_point() +
  labs(x = paste0("PC1 (",round(e$values[1],2),"%)"), 
       y = paste0("PC2 (",round(e$values[2],2),"%)"), 
       title = "TSHR region: 8.85-8.95Mb") +
  scale_color_manual(values = c("#3498db","#2ecc71","#f1c40f")) +
  theme_bw()
P2

#Write to file
write.table(select(PCs, sample_ID, cluster), "inferred_spawing_time.txt", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


###########
#Calculate HW

PCs$haplotype <- NA
for (ROW in 1:nrow(PCs)){
  if (!is.na(PCs$spawning.season.pc.based[ROW]) && PCs$spawning.season.pc.based[ROW] == "Spring-like"){
    PCs$haplotype[ROW] <- "0/0"
  } else if (!is.na(PCs$spawning.season.pc.based[ROW]) && PCs$spawning.season.pc.based[ROW] == "Intermediate"){
    PCs$haplotype[ROW] <- "0/1" 
  } else if (!is.na(PCs$spawning.season.pc.based[ROW]) && PCs$spawning.season.pc.based[ROW] == "Autumn-like"){
    PCs$haplotype[ROW] <- "1/1"
  } 
}

g1 <- genetics::genotype(as.vector(PCs$haplotype))
HWE <- genetics::HWE.test(g1)
HWE_statistic <- as.data.frame(HWE$test[3])




###########
#Plot female only phenos

#Load phenos
female.phenos.colnames <- colnames(read.table("../resources/curated_phenos_females.input.txt", header = TRUE))
female.phenos.data <- read.table("../resources/curated_phenos_females.input.txt", skip = 2, header = FALSE)
colnames(female.phenos.data) <- female.phenos.colnames

#Join phenos to cluster assignment
female.phenos.data <- female.phenos.data %>% dplyr::left_join(PCs)
female.phenos.data$cluster <- as.factor(female.phenos.data$cluster)

library(ggpubr)
compare_means(LC ~ cluster,  data = female.phenos.data)
my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )
ggplot(data = female.phenos.data, aes(x = cluster, y = LC, colour = cluster)) +
  geom_point(position = "jitter") +
  geom_boxplot(outlier.colour = NA, colour = "black", alpha = 0) +
  scale_color_manual(values = c("#3498db","#2ecc71","#f1c40f")) +
  stat_compare_means(comparisons = my_comparisons) +
  scale_x_discrete(labels = c("Spring-like","Intermediate","Autumn-like")) +
  theme(legend.position="none") +
  ggtitle("Leading cohort")



