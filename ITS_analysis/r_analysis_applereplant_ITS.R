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

library(BiocInstaller)
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
useDevel()
biocLite("microbiome")
library(microbiome)

# add rarefied otu table into phyloseq object
taxa_otu_table <- import_biom("/Users/arifinabintarti/Documents/Parent/apple_replant/ITS_analysis/single_rare.biom")
taxa_otu_table

#annotate otu table
rownames(map) <- map$sample_code
phyloseq_map <- sample_data(map)
taxa_otu_table_annotated <- merge_phyloseq(taxa_otu_table,phyloseq_map)
rank_names(taxa_otu_table_annotated)
colnames(tax_table(taxa_otu_table_annotated)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                                   o = "Order", f = "Family", g = "Genus", s = "Species")
taxa_otu_table_annotated

sample_phylum <- tax_glom(taxa_otu_table_annotated, taxrank = "Phylum", NArm = FALSE)
sample_phylumRA <- transform_sample_counts(sample_phylum, function(x)100* x/sum(x))
plot_bar(sample_phylumRA, fill="Phylum")                                                  
get_taxa_unique(sample_phylumRA, "Phylum")
sort(tapply(taxa_sums(sample_phylumRA), factor(tax_table(sample_phylumRA)[, "Phylum"]), sum), decreasing = TRUE)

#merge sample by site
site_phyl <- merge_samples(taxa_otu_table_annotated, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(taxa_otu_table_annotated)$site_name)
#merge taxa by rank
site_phylum <- tax_glom(site_phyl, taxrank = "Phylum", NArm = FALSE)
site_RA <- transform_sample_counts(site_phylum, function(x) 100 * x/sum(x))
#plot bar
plot_bar(site_RA, "site_name", fill="Phylum")

# Sum of Phylum
head(sort(tapply(taxa_sums(taxa_otu_table_annotated), factor(tax_table(taxa_otu_table_annotated)[, "Phylum"]), sum), decreasing = TRUE))


find.top.taxa <- function(x,taxa){
  require(phyloseq)
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa)) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- data.frame(tax[k,])
  m <- data.frame(otu[,k])
  s <- as.name(taxa)
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
find.top.taxa(top.pdy,"Phylum")



#samples name and variables
head(sample_names(taxa_otu_table_annotated))
#Total OTU abundance (reads) in each sample
sample_sums(taxa_otu_table_annotated)

meta <- meta(taxa_otu_table_annotated)
taxonomy <- tax_table(taxa_otu_table_annotated)
otu.absolute <- abundances(taxa_otu_table_annotated)
otu.relative <- abundances(taxa_otu_table_annotated, "compositional")
df <- psmelt(taxa_otu_table_annotated)
head(df)
topx <- top_taxa(taxa_otu_table_annotated, n = 10)
ranks <- rank_names(taxa_otu_table_annotated)
ranks
taxa  <- taxa(taxa_otu_table_annotated) 
taxa
head(get_taxa_unique(taxa_otu_table_annotated, "Phylum"))

###################calculate relative abundance along	taxonomical hierarchy################
fungi_phyla = tax_glom(taxa_otu_table_annotated, "Phylum", NArm = FALSE)
tax_table(fungi_phyla)
fungi_phyla_ab = transform_sample_counts(fungi_phyla, function(x) x/sum(x))
#or
abundances(fungi_phyla, "compositional")


##########calculate relative abundance along	taxonomical hierarchy in specific site#######
site_phyl <- merge_samples(taxa_otu_table_annotated, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(taxa_otu_table_annotated)$site_name)
site_fungi_phyla = tax_glom(site_phyl, "Phylum", NArm = FALSE)
site_fungi_phyla
site_fungi_phyla_ab = transform_sample_counts(site_fungi_phyla, function(x) x/sum(x))
tax_table(site_fungi_phyla_ab)
#or
abundances(site_fungi_phyla, "compositional")
tax_table(site_fungi_phyla)
otu_table(site_fungi_phyla)


taxa_otu_table_annotated
otu_table(taxa_otu_table_annotated)
sample_sums(taxa_otu_table_annotated)


otu_table(site_phyl)
otu.relative_each <- abundances(site_phyl, "compositional")
dim(otu.relative_each)

# 1. Clarksville
site_phyl
site_phyl_RA = transform_sample_counts(site_phyl, function(x) 100*x/sum(x))
site_phyl_RA
filter <- filter_taxa(site_phyl_RA, function(x) 100*x/sum(x), TRUE)
C <- sort(get_taxa(site_phyl_RA, "Clarksville"), decreasing = TRUE)[1:10]
C
# How to see the taxonomy 
C10 <- prune_taxa(names(C), site_phyl_RA)
C10
tax_table(C10)


