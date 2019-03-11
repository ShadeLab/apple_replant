#################################### Analysis ITS Sequences of 20 Sites only #####################################
# Date: May 30th 2018
# By : AF. Bintarti

# INSTALL PACKAGES
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
source("https://bioconductor.org/biocLite.R")
biocLite()
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/ITS_analysis/')
wd <- print(getwd())
otuITS <- read.table('single_rare45its.txt', sep='\t', header=T, row.names = 1)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
head(otuITS)
dim(otuITS)
taxonomy <- otuITS[,'taxonomy']
taxonomy
otuITS <- otuITS[,-46]
set.seed(13)
head(sort(colSums(otuITS, na.rm = FALSE, dims = 1), decreasing = FALSE))
sort(rowSums(otuITS, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(otuITS)
otuITS

# CALCULATE THE ALPHA DIVERSITY (SHANNON, RICHNESS, PIELOU)
otuITS_rare_PA <- 1*(otuITS>0)
sITS <- specnumber(otuITS, MARGIN = 2)
hITS <- diversity(t(otuITS), index = 'shannon')
pielouITS <- hITS/log(sITS)
map20 <- map[31:75,]
map.div <- map20
map.div$Richness <- sITS
map.div$Shannon <- hITS
map.div$Pielou <- pielouITS
names(map.div)
map.alphaITS <- melt(map.div, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness', 'Pielou'))

# GGPLOT OF THE ALPHA DIVERSITY (BASED ON 'site_name', 'cultivar', 'rootstock')
ggplot(map.alphaITS, aes(y=value, x=cultivar)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=cultivar))+
  geom_point(aes(color=cultivar))+
  theme(axis.text.x=element_text(angle=90, hjust=1))


ggplot(map.alphaITS, aes(y=value, x=site_name)) + 
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=site_name))+
  geom_point(aes(color=site_name))+
  theme(legend.position="bottom",axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=25),
        strip.text.x = element_text(size=30,colour = "black", face = "bold"),
        legend.text=element_text(size=20),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'),
        plot.title = element_text(size = rel(2)),axis.title=element_text(size=18,face="bold"))

# CALCULATE THE BETA DIVERSITY (PCA PLOT)
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_distITS <- vegdist(t(otuITS), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoaITS <- cmdscale(otu_distITS, eig=T)
env <- map20[,c(12:23, 25:37)]
# scores of PC1 and PC2
ax1ITS.scores=otu_pcoaITS$points[,1] 
ax2ITS.scores=otu_pcoaITS$points[,2] 
env_fitITS <- envfit(otu_pcoaITS, env, na.rm=TRUE)
ax1ITS <- otu_pcoaITS$eig[1]/sum(otu_pcoaITS$eig)
ax2ITS <- otu_pcoaITS$eig[2]/sum(otu_pcoaITS$eig)
map2=cbind(map20,ax1ITS.scores,ax2ITS.scores)
# PCoA plot with fit of environmental vectors
pcoa_plot <- plot(ax1ITS.scores, ax2ITS.scores, xlab=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))
plot(env_fitITS, p.max=0.05, col="red1")

###GGPLOT PCoA with Env_fit###
spp.scrsITS <- as.data.frame(scores(env_fitITS, display = "vectors"))
spp.scrsITS <- cbind(spp.scrsITS, env = rownames(spp.scrsITS))
#only significant pvalues
#shortcutting ef$vectors
A_ITS <- as.list(env_fitITS$vectors)
#creating the dataframe
pvals_ITS <- as.data.frame(A_ITS$pvals)
arrows_ITS <- as.data.frame(A_ITS$arrows*sqrt(A_ITS$r))
C_ITS <- cbind(arrows_ITS, pvals_ITS)
#subset
Cred_ITS<-subset(C_ITS,pvals_ITS<0.05)
Cred_ITS <- cbind(Cred_ITS, env = rownames(Cred_ITS))

ggplot(data = map2, aes(x=ax1ITS.scores, y=ax2ITS.scores, size=10)) + 
  geom_point(aes(color = site_name), size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  geom_segment(data=Cred_ITS,aes(x=0,xend=Dim1,y=0,yend=Dim2),
               arrow = arrow(length = unit(0.25, "cm")),size=0.5,colour="grey") + 
  geom_text(data=Cred_ITS,aes(x=Dim1,y=Dim2, label=env),size=5)+
  coord_fixed()+
  theme(legend.position="right",plot.title = element_text(size = rel(2), face="bold"),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=20),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))+
  ggtitle("ITS PCoA (Bray-Curtis)-Environmental fit")















# PCoA plot with symbol control in ggplot2
# plot by site_name
fig <- ggplot(data=map2, aes(x=ax1ITS.scores, y=ax2ITS.scores, size=8))+
  geom_point(aes(color=site_name)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig
# plot by rootstock
fig <- ggplot(data=map2, aes(x=ax1ITS.scores, y=ax2ITS.scores, size=8))+
  geom_point(aes(color=rootstock)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig
# plot by cultivar
fig <- ggplot(data=map2, aes(x=ax1ITS.scores, y=ax2ITS.scores, size=8))+
  geom_point(aes(color=cultivar)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig

# PERMANOVA
adonis(otu_distITS~map2$site_name)
adonis(otu_distITS~map2$rootstock)
adonis(otu_distITS~map2$cultivar)

# MAKE A BARPLOT USING PHYLOSEQ
# INSTALL PACKAGES
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)

# ADD THE TAXONOMY W/O PREFIX
# WRITE.CSV
write.csv(taxonomy, file = "its_taxonomy.csv")
# READ .CSV
tax_its = read.csv("its_taxon.csv", sep=',', header=T)
tax_its
# ADD OTU TABLE
dim(otuITS)
dim(tax_its)
rownames(tax_its) <- rownames(otuITS)
tax_its
class(otuITS)
class(tax_its)
OTU_its = otu_table(otuITS, taxa_are_rows = TRUE)
TAX_its = tax_table(as.matrix(tax_its))
# IMPORT OTU TABLE
otu_phyl <- import_biom("/Users/arifinabintarti/Documents/Parent/apple_replant/16S analysis/single_rare20.biom")
otu_phyl
# ADD MAP
rownames(map20) <- map20$sample_code
phyloseq_map <- sample_data(map20)
otu_map_its <- merge_phyloseq(OTU_its,TAX_its,phyloseq_map)
otu_map_its
# CHANGE THE RANK NAMES
rank_names(otu_map_its)
colnames(tax_table(otu_map)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
# MERGE SAMPLES BY SITE_NAMES
site_phyl <- merge_samples(otu_map_its, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(otu_map)$site_name)
# MERGE TAXA BY RANK
site_phylum <- tax_glom(site_phyl, taxrank = "Phylum", NArm=FALSE)
site_RA <- transform_sample_counts(site_phylum, function(x) 100 * x/sum(x))


# MAKE A BAR PLOT
plot_bar(site_RA, "site_name", fill="Phylum")+
  theme(axis.text.y=element_text(size=22),axis.text.x=element_text(size=22,angle=49,hjust =1, face = "bold"), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=16),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))


otu_table(site_RA)
otu_table(site_phylum)


phylum.sum = tapply(taxa_sums(site_phyl), tax_table(site_phyl)[, "Phylum"], sum, NArm=FALSE)

top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(site_phyl)[, "Phylum"] %in% top5phyla), site_phyl)
GP1
GP1_RA <- transform_sample_counts(GP1, function(x) 100 * x/sum(x))

otu_table(GP1_RA)
phylum.sum = tapply(taxa_sums(site_phylum), tax_table(site_phylum)[, "Phylum"], sum, NArm=FALSE)





library("RColorBrewer")
library("ggplot2")
library("plyr")
library("vegan")
library("reshape2")
library("ape")
library("phyloseq")
library("knitr")
library("xtable")
library("colorspace")
library("ggrepel")
#otu abundances and frequencies required for generating the basic plots
otu_map_its.ra = transform_sample_counts(otu_map_its, function(x) x/sum(x))
otu_table(otu_map_its.ra)
otu_map_its.ra

phylum.sum = tapply(taxa_sums(site_phyl), tax_table(site_phyl)[, "Phylum"], sum, NArm=FALSE)
top12phyla = names(sort(phylum.sum, TRUE))[1:12]

abund_val <- function(normalized){
  otu.abun = apply(otu_table(normalized),1,mean)
  otu.freq = rowSums(otu_table(normalized) != 0)/45
  phyla = as.vector(data.frame(tax_table(normalized))$Phylum)
  levels(phyla) = c(levels(phyla),"Other")
  keephyla = c("Ascomycota","Basidiomycota","Mortierellomycota", "Chytridiomycota", "Glomeromycota", "Rozellomycota", "Entomophthoromycota", "Mucoromycota", "Kickxellomycota", "Calcarisporiellomycota", "Blastocladiomycota")
  phyla[!(phyla %in% keephyla)] = "Other"
  phyla = as.vector(phyla)
  phyla=as.factor(phyla)
  otuabun = cbind.data.frame(abundance=log(otu.abun),frequency=otu.freq,phyla)
  return(otuabun)
}

abund_all <- abund_val(otu_map_its.ra)
# Use color brewer to pick a color scheme for the phyla
brew = brewer.pal(12, "Paired")
# Create a scatterplot of OTUs showing their average relative abundance and occupancy
ggplot(abund_all, aes(x=abundance,y=frequency,color=phyla)) + 
  geom_point(size=3) + xlab("Mean relative abundance (log10 scale)") + 
  ylab("Mean occupancy (n=45)") + scale_colour_brewer(palette="Paired")+
  labs(title="Fungal abundancy vs. occupancy plots")+ xlim(-14.6, -2)+
  theme(plot.title = element_text(size = rel(2)),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=18),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))


######## Occupancy-Abundancy #########

#Comulative Occ_Abund
otuITS_PA <- 1*((otuITS>0)==1)
otuITS_PA <- otuITS_PA[rowSums(otuITS_PA)>0,]
Occ_otuITS <- rowSums(otuITS_PA)/ncol(otuITS_PA)
abs_otuITS <- otuITS[rowSums(otuITS)>0,]
otuITS_rel <- decostand(abs_otuITS, method="total", MARGIN=2)
com_abund_otuITS <- rowSums(otuITS_rel)

#color code for the top most abundant OTUs
color_top <- com_abund_otuITS
color_top[] <- 'black' 
otuITS_top_abund <- log10(com_abund_otuITS)[order((com_abund_otuITS), decreasing=TRUE)][1:10] #top 10 most abundant OTUs 
color_top[names(color_top) %in% names(otuITS_top_abund)] <- 'red'
plot(log10(com_abund_otuITS), Occ_otuITS, col=color_top, pch=20, ylab='Occupancy (n=45)', xlab='log(sum relative abundace per OTU)\n (n=20838 OTUs)', 
     main='Fungal occupancy-abundance\n plot (comulative)')
#Mean Occ_Abund
otuITS_PA <- otuITS_PA[rowSums(otuITS_PA)>0,]
Mean_otuITS_PA <- apply(otuITS_PA, 1, mean)
Mean_Occ <- rowSums(otuITS_PA)/ncol(otuITS_PA)
Mean_abund <- apply(otuITS_rel, 1, mean)

#order(Mean_abund[names(Mean_abund) %in% unique(selected_otus_switch$otu)])

#Creating df for plotting with ggplot and adding color code for the shared and unique OTUs
df_occ <- data.frame(otu=names(Occ_otuITS), occ=Occ_otuITS) 
df_abun <- data.frame(otu=names(Mean_abund), abun=log10(Mean_abund))
otu_col <- data.frame(otu=names(color_top), col=color_top)
occ_abun <- left_join(df_abun, df_occ, by='otu')
occ_abun <- left_join(occ_abun, otu_col, by='otu')
occ_abun <- cbind.data.frame(occ_abun, phyla, by='otu')
occ_abun$unique <- 'unique'
#occ_abun$unique[occ_abun$otu %in% misc_occ_abun$otu] <- 'shared' #run after the miscanthus block bellow

#####################################
#Figure 3A - Occupancy abundance plot
#####################################
setEPS()
postscript('fungal_occ_abund.eps', width = 4.5, height = 4)
brew = brewer.pal(12, "Paired")
ggplot(data=occ_abun, aes(x=abun, y=occ, color=phyla)) +
  geom_point(size=4, pch=20)+
  scale_colour_brewer(palette="Paired")+
  labs(x="log(mean relative abundace per OTU)\n (n=2952 OTUs)", y= "Mean occupancy (n=45)", title='Fungal abundancy vs. occupancy plots') +
  geom_text_repel(data=occ_abun[occ_abun$col=='red',], aes(label=otu), box.padding = unit(0.45, "lines"), show.legend = FALSE) +
  geom_hline(aes(yintercept=.4), linetype='dashed', size=1.5) +
  geom_vline(aes(xintercept=-2.5), linetype='dashed', size=1.5) + 
  theme(plot.title = element_text(size = rel(2)),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"))+
  scale_y_continuous(breaks=seq(0,1,.2)) +
  scale_x_continuous(breaks=seq(-10,2,2))


dev.cur()
dev.off()



abund_all
phyla = as.vector(data.frame(tax_table(otu_map_its.ra))$Phylum)
phyla
levels(phyla) = c(levels(phyla),"Other")
keephyla = c("Ascomycota","Basidiomycota","Mortierellomycota", "Chytridiomycota", "Glomeromycota", "Rozellomycota", "Entomophthoromycota", "Mucoromycota", "Kickxellomycota", "Calcarisporiellomycota", "Blastocladiomycota")
phyla[!(phyla %in% keephyla)] = "Other"
phyla[!(phyla %in% keephyla)] = "Other"
phyla = as.vector(phyla)
phyla=as.factor(phyla)

occ_abun <- cbind.data.frame(occ_abun, phyla, by='otu')
occ_abun






