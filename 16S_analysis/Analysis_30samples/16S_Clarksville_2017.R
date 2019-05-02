##### Workflow for 16S sequence analysis data from 30 samples_Clarksville 2017 #####
# Date: March 28th 2019
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
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/16S_analysis/Analysis_30samples/')
wd <- print(getwd())
otu <- read.table('OTU_rarefied.txt', sep='\t', header=T, row.names = 1)
map <- read.table('AppleReplant_2017_Clarksville.csv', sep=',', header=TRUE)
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
otu <- otu[,-31]
class(otu)
names(otu)
dim(otu)
set.seed(201)
head(sort(colSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
dim(otu)
sum(otu)
###There are 1382/18655 Unassigned OTU for Clarksville 30 samples



# PHYLOSEQ
# INSTALL PACKAGES
source('http://bioconductor.org/biocLite.R')
install.packages("BiocManager")
library(BiocManager)
biocLite('phyloseq')
library(phyloseq)
# ADD THE TAXONOMY W/O PREFIX
# WRITE.CSV
write.csv(taxonomy, file = "16S_taxonomy.csv")
# READ .CSV
tax_16S = read.csv("16S_TAX.csv", sep=',', header=T)
tax_16S
# ADD OTU TABLE
dim(otu)
dim(tax_16S)
rownames(tax_16S) <- rownames(otu)
tax_16S
class(otu)
class(tax_16S)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_16S))
# ADD MAP
rownames(map) <- map$sample_code
phyloseq_map <- sample_data(map)
PHYL <- merge_phyloseq(OTU,TAX,phyloseq_map)
PHYL
# 1. Relative Abund across all 30 samples
# merge samples by "Phylum"
PHYLUM <- tax_glom(PHYL, taxrank = "Phylum", NArm=FALSE)
PHYLUM_RA <- transform_sample_counts(PHYLUM, function(x) 100 * x/sum(x))
# make barplot
plot_bar(PHYLUM_RA, fill="Phylum")+
 theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18,angle=49,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'))

TopNOTUs <- names(sort(taxa_sums(PHYLUM_RA), TRUE)[1:20])
ent20   <- prune_taxa(TopNOTUs, PHYLUM_RA)
plot_bar(ent20, fill="Phylum")+
 theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15,angle=90,hjust =1), axis.title=element_text(size=15),
       legend.text=element_text(size=10),legend.title = element_text(size = 10), legend.spacing.x = unit(0.1, 'cm'))
taxa_sums(ent20)/30

x <- sort(taxa_sums(ent20)/30, decreasing = T)
tax_table(ent20)

# RAREFACTION CURVE
rarecurve(t(otu_table(PHYL)), step=50, cex=0.5, xlab = "Reads", ylab = "Bacterial/Archaeal OTUs")

# CALCULATE THE BETA DIVERSITY (PCA PLOT)
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist <- vegdist(t(otu), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
otu_pcoa
env <- map[,c(12:23)]
# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1]
ax2.scores=otu_pcoa$points[,2] 
env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
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



















