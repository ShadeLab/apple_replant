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
otu_even26k.31_75 <- otu_even26k[,31:75]
otu_rare_PA <- 1*(otu_even26k>0)
otu_rare_PA.31_75 <- 1*(otu_even26k.31_75>0)
s <- specnumber(otu_even26k, MARGIN = 2)
s.31_75 <- specnumber(otu_even26k.31_75, MARGIN = 2)
h <- diversity(t(otu_even26k), index = 'shannon')
h.31_75 <- diversity(t(otu_even26k.31_75), index = 'shannon')
pielou <- h/log(s)
pielou.31_75 <- h.31_75/log(s.31_75)

map.20 <- map[31:75,]
map.div <- map
map.div.31_75 <- map[31:75,]
map.div$Richness <- s
map.div.31_75$Richness <- s.31_75
map.div$Shannon <- h
map.div.31_75$Shannon <- h.31_75
map.div$Pielou <- pielou
map.div.31_75$Pielou <- pielou.31_75
names(map.div)
names(map.div.31_75)
map.alpha <- melt(map.div, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness', 'Shannon', 'Pielou'))
map.alpha.20 <- melt(map.div.31_75, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness','Shannon', 'Pielou'))

ggplot(map.alpha, aes(y=value, x=site_name)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot()+
  geom_point()+
  theme(axis.text.x=element_text(angle=90, hjust=1))

# ggplot for 20 blocks
ggplot(map.alpha.20, aes(y=value, x=cultivar)) +
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=cultivar))+
  geom_point(aes(color=cultivar))+
  theme(axis.text.x=element_text(angle=90, hjust=1))

otu_dist <- vegdist(t(otu_even26k), method='bray')
otu_dist.20 <- vegdist(t(otu_even26k.31_75), method = 'bray')
otu_pcoa <- cmdscale(otu_dist, eig=T)
otu_pcoa.20 <- cmdscale(otu_dist.20, eig = T)
env=map[,c(9,11:22, 24:36)]
env.20=map[31:75,c(12:23, 25:37)]

ax1.scores=otu_pcoa$points[,1] 
ax2.scores=otu_pcoa$points[,2] 
ax1.scores.20=otu_pcoa.20$points[,1]
ax2.scores.20=otu_pcoa.20$points[,2]

env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
env_fit.20 <- envfit(otu_pcoa.20, env.20, na.rm = TRUE)
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax1.20 <- otu_pcoa.20$eig[1]/sum(otu_pcoa.20$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
ax2.20 <- otu_pcoa.20$eig[2]/sum(otu_pcoa.20$eig)
map2=cbind(map,ax1.scores,ax2.scores)
map.20 <- cbind(map.20, ax1.scores.20, ax2.scores.20)

#PCoA plot with fit of environmental vectors
pcoa_plot <- plot(otu_pcoa$points[,1], otu_pcoa$points[,2], xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fit, p.max=0.05, col="red")

pcoa_plot.20 <- plot(otu_pcoa.20$points[,1], otu_pcoa.20$points[,2], xlab=paste("PCoA1: ",round(ax1.20,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2.20,3)*100,"% var. explained", sep=""))
plot(env_fit.20, p.max=0.05, col="red")

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


###barplot###
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

#change the taxa's rank names
rank_names(taxa_otu_table)
colnames(tax_table(taxa_otu_table)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                  o = "Order", f = "Family", g = "Genus", s = "Species")
colnames(tax_table(taxa_otu_table_annotated_merge)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                                         o = "Order", f = "Family", g = "Genus", s = "Species")

#merge taxa by rank
taxa_otu_table_annotated_merge <- tax_glom(taxa_otu_table_annotated_merge, taxrank = "Phylum")

#plot bar

plot_bar(taxa_otu_table_annotated_merge, "site_name", fill="Phylum")

#barplot for top 10 Genus
TopNOTUs = names(sort(taxa_sums(taxa_otu_table_annotated_merge), decreasing=TRUE)[1:20])
print(TopNOTUs)
ent20 = prune_taxa(TopNOTUs, taxa_otu_table_annotated_merge)
ent20
plot_bar(ent20, "site_name", fill = "Phylum", facet_grid = ~Phylum) + 
  theme(axis.text.x=element_text(size=rel(0.5), angle=90))

p = plot_bar(ent10, "site_name", fill="Phylum", facet_grid=~Phylum)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

### Relative Abundance (Percentage) ###
colnames(tax_table(percent_abundance))
percent_abundance <- transform_sample_counts(ent20, function(x) 100 * x/sum(x))
percent_abundance_phylum <- tax_glom(percent_abundance, taxrank="Phylum")
percent_abundance.matrix = as(otu_table(percent_abundance_phylum), "matrix")

OTUdf = as.data.frame(percent_abundance.matrix)
OTUdf
genfac = factor(tax_table(ent20)[, "Phylum"])
genfac

gentab = apply(otu_table(ent20), MARGIN = 1, function(x) {
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

#for understanding which taxa are indicative of which rootstocks, etc
#indicspec (indicator species analysis)
#install indispec
install.packages("")
library(indicspecies)
indval_otu <- multipatt(t(otu_even26k), cluster=map$cultivar, control = how(nperm = 999))

##################################################################################################
###Phyloseq object for only 45 samples without Clarksville###

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

# add environmental map only for certain variables
phyloseq_map_var<- sample_data(map)[31:75,c(4,12:23, 25:37)]
phyloseq_map_var
taxa_otu_table1_20.taxamerged.var <- merge_phyloseq(taxa_otu_table1_20.taxamerged,phyloseq_map_var)
taxa_otu_table1_20.taxamerged.var


###PCoA Plot using Phyloseq###
install.packages("magrittr")
library(magrittr)
install.packages("ggplot2")
library(ggplot2)
packageVersion("vegan")
packageVersion("phyloseq")
packageVersion("ggplot2")
library(grid)
library(scales)
packageVersion("grid")

#Ordinate
taxa_pcoa <- ordinate(
  physeq = taxa_otu_table1_20.taxamerged.var, 
  method = "PCoA", 
  distance = "bray"
)
# Plot 
plot_ordination(
  physeq = taxa_otu_table1_20.taxamerged.var,
  ordination = taxa_pcoa,
  color = "site_name",
  title = "PCoA of Apple Orchard Bacterial Communities"
) +
  geom_point(aes(color = site_name), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

# Remove datapoints with NA
taxa_not_na <- taxa_otu_table1_20.taxamerged.var %>%
  subset_samples(
    !is.na(lime_index) 
  )
bray_not_na <- phyloseq::distance(physeq = taxa_not_na, method = "bray")
bray_not_na

# CAP ordinate
cap_ord <- ordinate(
  physeq = taxa_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ pH + lime_index + P_ppm + K_ppm + Ca_ppm + Mg_ppm + NO3N_ppm + NH4N_ppm 
  + OM_percent + sand_percent + silt_percent + clay_percent + Lesion + Dagger + Ring + Pin + Stunt + Spiral + Tylenchs + Aphelenchs + Dorylaims + Monochs + BacterialFeeders + MycorrhizalFungi + Oligochaetes
)
# CAP plot
cap_plot <- plot_ordination(
  physeq = taxa_not_na, 
  ordination = cap_ord, 
  color = "site_name", 
  axes = c(1,2)
) + 
  geom_point(aes(colour = site_name), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5)

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )



