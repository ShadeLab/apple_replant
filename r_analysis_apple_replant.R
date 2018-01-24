install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
source("https://bioconductor.org/biocLite.R")
biocLite()
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

setwd('/Users/arifinabintarti/Documents/apple_replant/')
otu <- read.table('Data/otu_table_v2.txt', sep='\t', header=T, row.names = 1)
map_16S <- read.csv('Data/map_16S.csv')
map_16S <- map_16S[-c(76,77),]
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
otu <- otu[,-76]

set.seed(13)
otu_rare <- t(rrarefy(t(otu), min(colSums(otu))))
head(otu_rare)

otu_rare <- otu_rare[rowSums(otu_rare)>0,]
otu_rare_PA <- 1*(otu_rare>0)
s <- specnumber(otu_rare, MARGIN = 2)
h <- diversity(t(otu_rare), index = 'shannon')
pielou <- h/log(s)

map.div <- map_16S
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou
names(map.div)
map.alpha <- melt(map.div, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness', 'Shannon', 'Pielou'))

ggplot(map.alpha, aes(y=value, x=site_name)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))

ggplot(map.alpha, aes(y=value, x=rootstock)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))


otu_dist <- vegdist(t(otu_rare), method='bray')
otu_pcoa <- cmdscale(otu_dist, eig=T)
env_fit <- envfit(otu_pcoa, map_16S[,c(11, 13, 14, 15, 16, 17, 18, 9, 4, 5, 6)])
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)

pcoa_plot <- plot(otu_pcoa$points[,1], otu_pcoa$points[,2], xlab='PCoA1 (14.7% var. explained)', ylab='PCoA2 (11.2% var. explained)')

adonis(otu_dist~map_16S$site_name)
  #Site location has a significant influence on composition
adonis(otu_dist~map_16S$rootstock)
adonis(otu_dist~map_16S$cultivar)
