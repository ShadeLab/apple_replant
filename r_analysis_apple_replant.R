#install.packages(c('vegan', 'tidyverse'))
#install.packages('reshape')
#source("https://bioconductor.org/biocLite.R")
#biocLite()
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

#setwd('/Users/arifinabintarti/Documents/apple_replant/')
otu <- read.table('otu_table_v2.txt', sep='\t', header=T, row.names = 1)
map <- read.table('clean_map_data.csv', header=TRUE)
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
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


#sample code for heatmap
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")
rowS=sort(rowSums(otu), decreasing=TRUE)
otu2=otu[names(rowS),]
dev.off()
#setEPS()
library(gplots)
#postscript("Figures/Fig5A.eps", width = 3.5, height=7, pointsize=10, paper="special")
heatmap.2(as.numeric(otu2[1:20,]),col=hc(100),scale="column",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", margins=c(5,13), srtCol=90)
#dev.off()

#for understanding which taxa are indicative of which rootstocks, etc
#indicspec (indicator species analysis)
#simper 
