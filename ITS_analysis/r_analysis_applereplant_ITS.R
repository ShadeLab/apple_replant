install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
source("https://bioconductor.org/biocLite.R")
biocLite()
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/ITS_analysis/')
wd <- print(getwd())
otu <- read.table('single_rare_its.txt', sep='\t', header=T, row.names = 1)
otu
colnames(otu)
dim(otu)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
taxonomy <- otu[,'taxonomy']
taxonomy
otu <- otu[,-76]
dim(otu)
set.seed(13)
head(sort(rowSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(colSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
otu_rare_PA <- 1*(otu>0)
s <- specnumber(otu, MARGIN = 2)
h <- diversity(t(otu), index = 'shannon')
pielou <- h/log(s)

map.div <- map
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou
names(map.div)
map.alpha <- melt(map.div, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness', 'Shannon', 'Pielou'))


ggplot(map.alpha, aes(y=value, x=rootstock)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=rootstock))+
  geom_point(aes(color=rootstock))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

otu_dist <- vegdist(t(otu), method='bray')
otu_pcoa <- cmdscale(otu_dist, eig=T)
env=map[,c(10,12:23, 25:37)]

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
  ggtitle("ITS PCoA (Bray-Curtis)")
fig

#plot by site_name
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=8))+
  geom_point(aes(color=site_name)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("ITS PCoA (Bray-Curtis)")
fig


adonis(otu_dist~map2$site_name)
adonis(otu_dist~map2$rootstock)
adonis(otu_dist~map2$cultivar)

########PHYLOSEQ################

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)

# add rarefied otu table into phyloseq object
taxa_otu_table <- import_biom("/Users/arifinabintarti/Documents/Parent/apple_replant/ITS analysis/single_rare.biom")
taxa_otu_table

#annotate otu table
rownames(map) <- map$sample_code
phyloseq_map <- sample_data(map)
taxa_otu_table_annotated <- merge_phyloseq(taxa_otu_table,phyloseq_map)
rank_names(taxa_otu_table_annotated)
colnames(tax_table(taxa_otu_table_annotated)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                                   o = "Order", f = "Family", g = "Genus", s = "Species")


sample_phylum <- tax_glom(taxa_otu_table_annotated, taxrank = "Phylum")
sample_phylumRA <- transform_sample_counts(sample_phylum, function(x)100* x/sum(x))
plot_bar(sample_phylumRA, fill="Phylum")                                                  
get_taxa_unique(sample_phylumRA, "Phylum")
sort(tapply(taxa_sums(sample_phylumRA), factor(tax_table(sample_phylumRA)[, "Phylum"]), sum), decreasing = TRUE)

#merge sample by site
site_phyl <- merge_samples(taxa_otu_table_annotated, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(taxa_otu_table_annotated)$site_name)
#merge taxa by rank
site_phylum <- tax_glom(site_phyl, taxrank = "Phylum")
site_RA <- transform_sample_counts(site_phylum, function(x) 100 * x/sum(x))
#plot bar
plot_bar(site_RA, "site_name", fill="Phylum")

# Sum of Phylum
head(sort(tapply(taxa_sums(taxa_otu_table_annotated), factor(tax_table(taxa_otu_table_annotated)[, "Phylum"]), sum), decreasing = TRUE))
