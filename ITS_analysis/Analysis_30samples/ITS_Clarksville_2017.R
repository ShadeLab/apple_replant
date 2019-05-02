##### Workflow for ITS sequence analysis data from 30 samples_Clarksville 2017 #####
# Date: March 30th 2019
# By : A. Fina Bintarti

# INSTALL PACKAGES
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
source("https://bioconductor.org/biocLite.R")
biocLite()
install.packages('BiocManager')
install.packages("ggrepel")
library(ggrepel)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)
library(devtools)
install.packages("agricolae")
library(agricolae)

#SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/ITS_analysis/Analysis_30samples/')
wd <- print(getwd())
otuITS <- read.table('OTU_ITS_rarefied.txt', sep='\t', header=T, row.names = 1)
map <- read.table('AppleReplant_2017_Clarksville.csv', sep=',', header=TRUE)
head(otuITS)
dim(otuITS)
otuITS_PA <- 1*(otuITS>0)
set.seed(201)
head(sort(colSums (otuITS, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums (otuITS, na.rm = FALSE, dims = 1), decreasing = FALSE))

# PHYLOSEQ
# INSTALL PACKAGES
source('http://bioconductor.org/biocLite.R')
install.packages("BiocManager")
library(BiocManager)
biocLite('phyloseq')
library(phyloseq)
# READ OTU_TABLE
otu_table = read.csv("ITS_OTU_table_Clarksville.csv", sep=",", row.names=1)
otu_table = as.matrix(otu_table)
dim(otu_table)
head(sort(colSums(otu_table, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otu_table, na.rm = FALSE, dims = 1), decreasing = FALSE))
# READ TAXONOMY
# seperated by kingdom phylum class order family genus species 
tax = read.csv("ITS_taxonomy_Clarksville.csv", sep=",", row.names=1)
tax = as.matrix(tax)
dim(tax)
tax
# IMPORT "otu_table", "tax", AND "map20" FILES INTO PHYLOSEQ OBJECT
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)
OTU
TAX <-  tax_table(tax)
TAX
rownames(map) <- map$sample_code
MAP <- sample_data(map)
# merge into one phyloseq object
physeqITS <- phyloseq(OTU,TAX,MAP)
physeqITS.ra <- transform_sample_counts(physeqITS, function(x) x/sum(x))

# 1. Relative Abund across all 30 samples
# merge samples by "Phylum"
PHYLUM <- tax_glom(physeqITS, taxrank = "Phylum", NArm=FALSE)
PHYLUM_RA <- transform_sample_counts(PHYLUM, function(x) 100 * x/sum(x))
# merge samples by "Class"
CL <- tax_glom(physeqITS, taxrank = "Class", NArm=FALSE)
CL_RA <- transform_sample_counts(CL, function(x) 100 * x/sum(x))
plot_bar(CL_RA, fill="Class")+
 theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18,angle=90,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'))
TopNOTUs <- names(sort(taxa_sums(CL_RA), TRUE)[1:20])
ent20   <- prune_taxa(TopNOTUs, CL_RA)
plot_bar(ent20, fill="Class")+theme_bw()+
 theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15,angle=90,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# make barplot
plot_bar(PHYLUM_RA, fill="Phylum")+
 theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18,angle=90,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'))

TopNOTUs <- names(sort(taxa_sums(PHYLUM_RA), TRUE)[1:10])
ent10   <- prune_taxa(TopNOTUs, PHYLUM_RA)
plot_bar(ent10, fill="Phylum")+theme_bw()+
 theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15,angle=90,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

taxa_sums(ent10)/30

x <- sort(taxa_sums(ent10)/30, decreasing = T)
tax_table(ent10)

# MERGE SAMPLES BY Sample Location
site_phyl <- merge_samples(physeqITS, "Sample.Location")
# merge samples by "Phylum"
site_PHYLUM <- tax_glom(site_phyl, taxrank = "Phylum", NArm=FALSE)
site_PHYLUM_RA <- transform_sample_counts(site_PHYLUM, function(x) 100 * x/sum(x))
plot_bar(site_PHYLUM_RA, fill="Phylum")+
 theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18,angle=90,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'))

DR <- subset_samples(physeqITS, sample_code%in%c("C25", "C26", "C27", "C28", "C29","C30"))
DR_phyl <- tax_glom(DR, taxrank = "Phylum", NArm=FALSE)
DR_ra <- transform_sample_counts(DR_phyl, function(x) 100 * x/sum(x))
TopNOTUs <- names(sort(taxa_sums(DR_ra), TRUE)[1:10])
ent10   <- prune_taxa(TopNOTUs, DR_ra)
taxa_sums(ent10)/6
x <- sort(taxa_sums(ent10)/6, decreasing = T)
tax_table(ent10)

RZ <- subset_samples(physeqITS, sample_code%in%c("C01", "C02", "C03", "C04", "C05","C06","C07","C08","C09","C10",
                                                 "C11","C12","C13","C14","C15","C16","C17","C18","C19","C20",
                                                 "C21","C22","C23","C24"))
RZ_phyl <- tax_glom(RZ, taxrank = "Phylum", NArm=FALSE)
RZ_ra <- transform_sample_counts(RZ_phyl, function(x) 100 * x/sum(x))
TopNOTUs <- names(sort(taxa_sums(RZ_ra), TRUE)[1:10])
ent10   <- prune_taxa(TopNOTUs, RZ_ra)
taxa_sums(ent10)/24
x <- sort(taxa_sums(ent10)/24, decreasing = T)
tax_table(ent10)







# RAREFACTION CURVE
rarecurve(t(otu_table(physeqITS)), step=50, cex=0.5, xlab = "Reads", ylab = "Fungal OTUs")

# CALCULATE THE BETA DIVERSITY (PCA PLOT)
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_distITS <- vegdist(t(otuITS), method='bray')
adonis(otu_distITS~map2$Sample.Location)
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoaITS <- cmdscale(otu_distITS, eig=T)
otu_pcoaITS
env <- map[,c(12:23)]
# scores of PC1 and PC2
ax1.scores=otu_pcoaITS$points[,1]
ax2.scores=otu_pcoaITS$points[,2] 
env_fit <- envfit(otu_pcoaITS, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1 <- otu_pcoaITS$eig[1]/sum(otu_pcoaITS$eig)
ax2 <- otu_pcoaITS$eig[2]/sum(otu_pcoaITS$eig)
map2=cbind(map,ax1.scores,ax2.scores)

pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fit, p.max=0.05, col="red1")

#shortcutting ef$vectors
A <- as.list(env_fit$vectors)
#creating the dataframe
pvals<-as.data.frame(A$pvals)
# environment scores (vectors scaled by R2 values)
clark.scores1 < -as.data.frame(scores(env_fit, display="vectors"))
clark.scores2 <- cbind(clark.scores1, pvals)
clark.scores3 <- cbind(clark.scores2,Variable=rownames(clark.scores2))
clark.scores4 <- subset(clark.scores3,pvals<0.05)

mult <-.25
(Fig1<-ggplot(data = map2, aes(x=ax1.scores, y=ax2.scores)) +
  theme_bw()+
  geom_point(aes(ax1.scores, y=ax2.scores,color=Sample.Location),size=2.5,shape=20,stroke=1.75)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  coord_fixed()+
  geom_segment(data=clark.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = clark.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
  theme(text = element_text(size=12),legend.title=element_text(size=12), 
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

#ggsave("Figure1.pdf",Fig1,width=190,height=120,units="mm")
ggsave("Figure1.tiff",Fig1,width=6, height=4, units="in",dpi=300)


