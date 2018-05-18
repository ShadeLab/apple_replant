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

#setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/16S analysis')
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/16S analysis')
#print working directory for future references
wd <- print(getwd())
otu <- read.table('single_rare1.txt', sep='\t', header=T, row.names = 1)
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

head(sort(colSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
#otu_even26k <- t(rrarefy(t(otu), min(colSums(otu))))
#head(otu_even26k)
#head(sort(rowSums (otu_even26k, na.rm = FALSE, dims = 1), decreasing = FALSE))
#summary(otu_even26k)
sum(otu)
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


ggplot(map.alpha, aes(y=value, x=cultivar)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=cultivar))+
  geom_point(aes(color=cultivar))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

#dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist <- vegdist(t(otu), method='bray')
#CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
env=map[,c(10,12:23, 25:37)]
#scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1] 
ax2.scores=otu_pcoa$points[,2] 

env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
map2=cbind(map,ax1.scores,ax2.scores)

#PCoA plot with fit of environmental vectors

pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fit, p.max=0.05, col="red1")


#PCoA plot with symbol control in ggplot2
#colors=c( "darkolivegreen1", "darkolivegreen4","burlywood", "burlywood4")
#plot by site_name
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=8))+
  geom_point(aes(color=site_name)) + 
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
#PERMANOVA
adonis(otu_dist~map2$site_name)
adonis(otu_dist~map2$rootstock)
adonis(otu_dist~map2$cultivar)

############################################################  
##########heatmap###############
install.packages("gplots")
library(gplots)
heatmap(as.matrix(otu),col=hc(100),scale="column",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", margins=c(5,13), srtCol=90)
###########################################################

#barplot using phyloseq
install.packages("phyloseq")
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)

# add rarefied otu table into phyloseq object
taxa_otu_table <- import_biom("/Users/arifinabintarti/Documents/Parent/apple_replant/16S analysis/single_rare.biom")
taxa_otu_table

#annotate otu table
rownames(map) <- map$sample_code
phyloseq_map <- sample_data(map)
taxa_otu_table_annotated <- merge_phyloseq(taxa_otu_table,phyloseq_map)
rank_names(taxa_otu_table_annotated)
colnames(tax_table(taxa_otu_table_annotated)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
 

sample_phylum <- tax_glom(taxa_otu_table_annotated, taxrank = "Phylum")
sample_phylumRA <- transform_sample_counts(sample_phylum, function(x)100* x/sum(x))
plot_bar(sample_phylumRA, fill="Phylum")                                                  
get_taxa_unique(sample_phylumRA, "Phylum")
sort(tapply(taxa_sums(sample_phylumRA), factor(tax_table(sample_phylumRA)[, "Phylum"]), sum), decreasing = TRUE)

##############################################################################                                                                                                          o = "Order", f = "Family", g = "Genus", s = "Species")
#merge sample by site
site_phyl <- merge_samples(taxa_otu_table_annotated, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(taxa_otu_table_annotated)$site_name)
#merge taxa by rank
site_phylum <- tax_glom(site_phyl, taxrank = "Phylum")
site_RA <- transform_sample_counts(site_phylum, function(x) 100 * x/sum(x))
#plot bar
plot_bar(site_RA, "site_name", fill="Phylum")
###############################################################################
#top20Phyla
site_phylum <- tax_glom(site_phyl, taxrank = "Phylum", NArm = FALSE)
TopNOTUs <- names(sort(taxa_sums(site_phylum), TRUE)[1:20])
site_RA <- transform_sample_counts(site_phylum, function(x) 100 * x/sum(x))
x4 <- prune_taxa(TopNOTUs, site_RA)
plot_bar(x4, fill="Phylum")
#############################################################################

phyloGlom = tax_glom(taxa_otu_table_annotated, "Phylum")
glomTax = tax_table(phyloGlom)[,"Phylum"]
glomOTU = otu_table(phyloGlom)
glomTable = merge(glomOTU,glomTax,by=0,all=TRUE)
rownames(glomTable) = glomTable[,"Phylum"]
glomTable$Row.names = NULL
glomTable$Phylum = NULL
head(glomTable)

colSums(glomTable)
rowSums(glomTable)
nrow(glomTable)

glomTable2 = glomTable / rep(colSums(glomTable), each = nrow(glomTable))
glomTable2

sort(rowSums(glomTable), decreasing = TRUE)
sort(rowSums(glomTable2), decreasing = TRUE)

##############################################################################
### barplot for top 10 Phylum ###
TopNOTUs = names(sort(taxa_sums(taxa_otu_table_annotated), decreasing=TRUE)[1:10])
print(TopNOTUs)
ent10 = prune_taxa(TopNOTUs, taxa_otu_table_annotated)
ent10
plot_bar(ent10, "site_name", fill = "Phylum", facet_grid = ~Phylum) + 
  theme(axis.text.x=element_text(size=rel(0.5), angle=90))

p = plot_bar(ent10, "site_name", fill="Phylum", facet_grid=~Phylum)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

### Relative Abundance (Percentage) ###
colnames(tax_table(percent_abundance))
percent_abundance <- transform_sample_counts(ent10, function(x) 100 * x/sum(x))
percent_abundance_phylum <- tax_glom(percent_abundance, taxrank="Phylum")
percent_abundance.matrix = as(otu_table(percent_abundance_phylum), "matrix")

OTUdf = as.data.frame(percent_abundance.matrix)
OTUdf
genfac = factor(tax_table(ent10)[, "Phylum"])
genfac

gentab = apply(otu_table(ent10), MARGIN = 1, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
gentab


#boxplot environmental variables

# 1. pH
pH.means <- aggregate(pH ~  site_name, map.div.31_75, mean)
ggplot(data=map.31_75, aes(x=site_name, y=pH, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = pH.means, aes(label = round(pH, digits = 2), y = pH + 0.15))

# 1. pH
pH.means <- aggregate(pH ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=pH, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = pH.means, aes(label = round(pH, digits = 2), y = pH + 0.15))
# 2. P_ppm
P_ppm.means <- aggregate(P_ppm ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=P_ppm, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = P_ppm.means, aes(label = round(P_ppm, digits = 2), y = P_ppm + 12.7))
# 3. K_ppm
K_ppm.means <- aggregate(K_ppm ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=K_ppm, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = K_ppm.means, aes(label = round(K_ppm, digits = 2), y = K_ppm + 17))
# 4. Ca_ppm
Ca_ppm.means <- aggregate(Ca_ppm ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=Ca_ppm, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = Ca_ppm.means, aes(label = round(Ca_ppm, digits = 2), y = Ca_ppm + 100))
# 5. Mg_ppm
Mg_ppm.means <- aggregate(Mg_ppm ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=Mg_ppm, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = Mg_ppm.means, aes(label = round(Mg_ppm, digits = 2), y = Mg_ppm + 15))
# 6. NO3N_ppm
NO3N_ppm.means <- aggregate(NO3N_ppm ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=NO3N_ppm, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = NO3N_ppm.means, aes(label = round(NO3N_ppm, digits = 2), y = NO3N_ppm + 0.1))
# 7. NH4N_ppm
NH4N_ppm.means <- aggregate(NH4N_ppm ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=NH4N_ppm, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = NH4N_ppm.means, aes(label = round(NH4N_ppm, digits = 2), y = NH4N_ppm + 0.48))
# 8. OM_percent
OM_percent.means <- aggregate(OM_percent ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=OM_percent, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = OM_percent.means, aes(label = round(OM_percent, digits = 2), y = OM_percent + 0.18))
# 9. sand_percent
sand_percent.means <- aggregate(sand_percent ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=sand_percent, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = sand_percent.means, aes(label = round(sand_percent, digits = 2), y = sand_percent + 2.5))
# 10. silt_percent
silt_percent.means <- aggregate(silt_percent ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=silt_percent, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = silt_percent.means, aes(label = round(silt_percent, digits = 2), y = silt_percent + 2.5))
# 11. clay_percent
clay_percent.means <- aggregate(clay_percent ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=clay_percent, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = clay_percent.means, aes(label = round(clay_percent, digits = 2), y = clay_percent + 1))
# 11. clay_percent
clay_percent.means <- aggregate(clay_percent ~  site_name, map.div.31_75, mean)
ggplot(data=map.div.31_75, aes(x=site_name, y=clay_percent, fill=site_name)) + geom_boxplot() + 
  scale_x_discrete(label=abbreviate) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = clay_percent.means, aes(label = round(clay_percent, digits = 2), y = clay_percent + 1))

##################################################################################################
#for understanding which taxa are indicative of which rootstocks, etc
#indicspec (indicator species analysis)
#install indispec
install.packages("")
library(indicspecies)
indval_otu <- multipatt(t(otu_even26k), cluster=map$cultivar, control = how(nperm = 999))

##################################################################################################
###Phyloseq object for only 45 samples without Clarksville###
#################################################################################################
taxa_otu_table
taxa_otu_table1_20 <- otu_table(taxa_otu_table)[,31:75]
taxa_otu_table1_20.taxamerged <- merge_phyloseq(taxa_otu_table1_20,tax_table(taxa_otu_table))
taxa_otu_table1_20.taxamerged
# add environmental map
rownames(map) <- map$sample_code
phyloseq_map1_20 <- sample_data(map)[31:75,]
phyloseq_map1_20
taxa_otu_table1_20.taxamerged.phyl <- merge_phyloseq(taxa_otu_table1_20.taxamerged,phyloseq_map1_20)
taxa_otu_table1_20.taxamerged.phyl

# relative abundances per 
abundance_phyla = tax_glom(taxa_otu_table1_20.taxamerged.phyl, "Phylum")
abundance_phyla
abundance_phyla_NArmFALSE = tax_glom(taxa_otu_table1_20.taxamerged.phyl, "Phylum", NArm = FALSE)
abundance_phyla_NArmFALSE
ntaxa(taxa_otu_table1_20.taxamerged.phyl); ntaxa(abundance_phyla)
abundance_genus = tax_glom(taxa_otu_table1_20.taxamerged.phyl, "Genus")
abundance_genus
ntaxa(taxa_otu_table1_20.taxamerged.phyl); ntaxa(abundance_genus)

tax_table(_phyla)






