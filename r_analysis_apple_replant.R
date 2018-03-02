#install.packages(c('vegan', 'tidyverse'))
install.packages(c('vegan', 'tidyverse'))
#install.packages('reshape')
install.packages('reshape')
#source("https://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
#biocLite()
biocLite()
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

#setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/')
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/')
#print working directory for future references
wd <- print(getwd())
otu <- read.table('otu_table_v2.txt', sep='\t', header=T, row.names = 1)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
taxonomy
otu <- otu[,-76]

set.seed(13)
#summary(otu)
#min. number of sequences is ~26K
#subsample to an even number of sequences per sample

otu_even26k <- t(rrarefy(t(otu), min(colSums(otu))))
head(otu_even26k)
#summary(otu_even26k)

otu_even26k <- otu_even26k[rowSums(otu_even26k)>0,]
otu_rare_PA <- 1*(otu_even26k>0)
s <- specnumber(otu_even26k, MARGIN = 2)
h <- diversity(t(otu_even26k), index = 'shannon')
pielou <- h/log(s)

map.div <- map
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou
names(map.div)
map.alpha <- melt(map.div, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness', 'Shannon', 'Pielou'))

ggplot(map.alpha, aes(y=value, x=site_name)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot()+
  geom_point()+
  theme(axis.text.x=element_text(angle=90, hjust=1))

ggplot(map.alpha, aes(y=value, x=rootstock)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot()+
  geom_point()+
  theme(axis.text.x=element_text(angle=90, hjust=1))

ggplot(map.alpha, aes(y=value, x=cultivar)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot()+
  geom_point()+
  theme(axis.text.x=element_text(angle=90, hjust=1))

otu_dist <- vegdist(t(otu_even26k), method='bray')
otu_pcoa <- cmdscale(otu_dist, eig=T)
env=map[,c(9,11:22, 24:36)]
ax1.scores=otu_pcoa$points[,1] 
ax2.scores=otu_pcoa$points[,2] 

env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
map2=cbind(map,ax1.scores,ax2.scores)

#PCoA plot with fit of environmental vectors
pcoa_plot <- plot(otu_pcoa$points[,1], otu_pcoa$points[,2], xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fit, p.max=0.05, col="red")

#PCoA plot with symbol control in ggplot2
#colors=c( "darkolivegreen1", "darkolivegreen4","burlywood", "burlywood4")
#plot by rootstock
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=8))+
  geom_point(aes(color=rootstock)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig

#plot by cultivar
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=8))+
  geom_point(aes(color=cultivar)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig


adonis(otu_dist~map2$site_name)
adonis(otu_dist~map2$rootstock)
adonis(otu_dist~map2$cultivar)
  
#heatmap
install.packages("gplots")
library(gplots)

#select most abundant OTUs
hc <- colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")
rowS_ordered <- sort(rowSums(otu_even26k), decreasing=TRUE)
otu2 <- otu_even26k[names(rowS_ordered),]
dev.off()

#dev.off()
#setEPS()
#postscript("Figures/Fig5A.eps", width = 3.5, height=7, pointsize=10, paper="special")
heatmap.2(as.matrix(otu2[1:20,]),col=hc(100),scale="column",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", margins=c(5,13), srtCol=90)


#barplot
install.packages("phyloseq")
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
taxa_otu_table <- import_biom("/Users/arifinabintarti/Documents/Parent/Data/single_rare.biom")
taxa_otu_table

#annotate otu table
rownames(map) <- map$sample_code
phyloseq_map <- sample_data(map)
taxa_otu_table_annotated <- merge_phyloseq(taxa_otu_table,phyloseq_map)

#merge sample by site
taxa_otu_table_annotated_merge <- merge_samples(taxa_otu_table_annotated, "site_name")
sample_data(taxa_otu_table_annotated_merge)$site_name <- levels(sample_data(taxa_otu_table_annotated)$site_name)
taxa_otu_table_annotated_merge <- transform_sample_counts(taxa_otu_table_annotated_merge, function(x) 100 * x/sum(x))

#merge sample by cultivar
taxa_otu_table_annotated_merge <- merge_samples(taxa_otu_table_annotated, "cultivar")
sample_data(taxa_otu_table_annotated_merge)$cultivar <- levels(sample_data(taxa_otu_table_annotated)$cultivar)
taxa_otu_table_annotated_merge <- transform_sample_counts(taxa_otu_table_annotated_merge, function(x) 100 * x/sum(x))

#merge sample by rootstock
taxa_otu_table_annotated_merge <- merge_samples(taxa_otu_table_annotated, "rootstock")
sample_data(taxa_otu_table_annotated_merge)$rootstock <- levels(sample_data(taxa_otu_table_annotated)$rootstock)
taxa_otu_table_annotated_merge <- transform_sample_counts(taxa_otu_table_annotated_merge, function(x) 100 * x/sum(x))

#merge sample by sample_location
taxa_otu_table_annotated_merge <- merge_samples(taxa_otu_table_annotated, "sample_location")
sample_data(taxa_otu_table_annotated_merge)$sample_location <- levels(sample_data(taxa_otu_table_annotated)$sample_location)
taxa_otu_table_annotated_merge <- transform_sample_counts(taxa_otu_table_annotated_merge, function(x) 100 * x/sum(x))
#merge taxa by rank
taxa_otu_table_annotated_merge <- tax_glom(taxa_otu_table_annotated_merge, taxrank = "Rank2")

#plot bar
plot_bar(taxa_otu_table_annotated_merge, "site_name", fill="Rank2")
plot_bar(taxa_otu_table_annotated_merge, "cultivar", fill="Rank2")
plot_bar(taxa_otu_table_annotated_merge, "rootstock", fill="Rank2")
plot_bar(taxa_otu_table_annotated_merge, "sample_location", fill="Rank2")


#for understanding which taxa are indicative of which rootstocks, etc
#indicspec (indicator species analysis)
 
#install indispec
install.packages("")
library(indicspecies)
indval_otu <- multipatt(t(otu_even26k), cluster=map$cultivar, control = how(nperm = 999))













