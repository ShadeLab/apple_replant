#################################### Analysis 16S Sequences of 20 Sites only #####################################
# Date: May 29th 2018
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
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/16S analysis/Analysis_45samples')
wd <- print(getwd())
otu <- read.table('otu_table_20.txt', sep='\t', header=T, row.names = 1)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
taxonomy
otu <- otu[,-46]
set.seed(13)
head(sort(colSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
dim(otu)
sum(otu)

# CALCULATE THE ALPHA DIVERSITY (SHANNON, RICHNESS, PIELOU)
otu_rare_PA <- 1*(otu>0)
s <- specnumber(otu, MARGIN = 2)
h <- diversity(t(otu), index = 'shannon')
pielou <- h/log(s)
map20 <- map[31:75,]
map.div <- map20
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou
names(map.div)
map.alpha <- melt(map.div, id.vars=c('site_name', 'cultivar', 'rootstock'), measure.vars=c('Richness', 'Pielou'))

# GGPLOT OF THE ALPHA DIVERSITY (BASED ON 'site_name', 'cultivar', 'rootstock')
ggplot(map.alpha, aes(y=value, x=site_name)) + 
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=site_name))+
  geom_point(aes(color=site_name))+
  theme(axis.text.x = element_text(size=12,angle=90,hjust =1),strip.text.x = element_text(size=22,colour = "black", face = "bold"),
  plot.title = element_text(size = rel(2)),axis.title=element_text(size=18,face="bold"))
      
ggplot(map.alpha, aes(y=value, x=site_name)) + 
  facet_wrap(~variable, scales = 'free_y') +
  geom_boxplot(aes(color=site_name))+
  geom_point(aes(color=site_name))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=25),
        strip.text.x = element_text(size=30,colour = "black", face = "bold"),
        legend.text=element_text(size=20),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'),
        plot.title = element_text(size = rel(2)),axis.title=element_text(size=18,face="bold"))







                                  

# CALCULATE THE BETA DIVERSITY (PCA PLOT)
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist <- vegdist(t(otu), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
env <- map20[,c(12:23, 25:37)]
# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1] 
ax2.scores=otu_pcoa$points[,2] 
env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
map2=cbind(map20,ax1.scores,ax2.scores)

###GGPLOT PCoA with Env_fit###
spp.scrs <- as.data.frame(scores(env_fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, env = rownames(spp.scrs))
#only significant pvalues
#shortcutting ef$vectors
A <- as.list(env_fit$vectors)
#creating the dataframe
pvals<-as.data.frame(A$pvals)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
C<-cbind(arrows, pvals)
#subset
Cred<-subset(C,pvals<0.05)
Cred <- cbind(Cred, env = rownames(Cred))

ggplot(data = map2, aes(x=ax1.scores, y=ax2.scores, size=10)) + 
  geom_point(aes(color = site_name), size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  geom_segment(data=Cred,aes(x=0,xend=Dim1,y=0,yend=Dim2),
               arrow = arrow(length = unit(0.25, "cm")),size=0.5,colour="grey") + 
  geom_text(data=Cred,aes(x=Dim1,y=Dim2, label=env),size=5)+
  coord_fixed()+
  theme(plot.title = element_text(size = rel(2), face="bold"),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
                     legend.text=element_text(size=20),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))+
  ggtitle("16S rRNA PCoA (Bray-Curtis)-Environmental fit")

# PCoA plot with fit of environmental vectors
pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fit, p.max=0.05, col="red1")
# add text into the plot
plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""), 
     main="16S rRNA PCoA (Bray-Curtis)",	type="n")
text(ax1.scores, ax2.scores, labels = row.names(t(otu)), cex=.7)
plot(env_fit, p.max=0.05, col="red1")



##################################################################
map20
unique(map20$site_name)
site=rep('black',nrow(map20))
site[map20$site_name=="Vinton(Schwallier)"]='deepskyblue'
site[map20$site_name=="SchweitzerBlock1,Zone15"]='deeppink'
site[map20$site_name=="SchweitzerBlock2,Zone3"]='darkblue'
site[map20$site_name=="Thome"]='cyan'
site[map20$site_name=="Dunnabeck"]='darkslategray3'
site[map20$site_name=="KlenkWest"]='magenta'
site[map20$site_name=="KlenkEast"]='darkorchid3'
site[map20$site_name=="DietrichWest"]='red'
site[map20$site_name=="DietrichEast "]='#a65628'
site[map20$site_name=="Wells"]='darkgoldenrod'
site[map20$site_name=="FrensNorth"]='4daf4a'
site[map20$site_name=="FrensSouth"]='#1919ff'
site[map20$site_name=="RoosinckPatterson"]='#ff7f00'
site[map20$site_name=="RoosinckReeman"]='#ffff33'
site[map20$site_name=="Rasch"]='coral1'
site[map20$site_name=="RothNorth"]='#f781bf'
site[map20$site_name=="RothSouth"]='green'
site[map20$site_name=="Kropf"]='blue4'
site[map20$site_name=="RaschBigRock"]='cornflowerblue'
site[map20$site_name=="RaschSouthofHouse"]='#984ea3'
with(map20,levels(site_name))
scl <- 3
plot(ax1.scores, ax2.scores,main="16S rRNA PCoA (Bray-Curtis)", xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
with(map2,points(ax1.scores, ax2.scores, display="site_name",col = values[site],
                 scaling=scl, pch=21, bg=values[site]))
text(ax1.scores, ax2.scores,display="sample_code", scaling=scl,
     cex=0.8, col = "black")
with(map2, legend("topright", legend=levels(site_name), bty="n",
                  col=values, pch=21, pt.bg=values))
legend('bottomright',c('Vinton(Schwallier)','SchweitzerBlock1,Zone15','SchweitzerBlock2,Zone','Thome','Dunnabeck','KlenkWest','KlenkEast','DietrichWest','DietrichEast','Wells','FrensNorth','FrensSouth','RoosinckPatterson','RoosinckReeman','Rasch','RothNorth','RothSouth','Kropf','RaschBigRock','RaschSouthofHouse'),
       pch=21,pt.bg=colvec,lty=0)
plot(env_fit, p.max=0.10, col="black", cex=1)

with(map20, levels(site_name))
library(randomcoloR)
n <- 20
values = c("#a65628", "red", "darkslategray3",
           "#4daf4a", "#1919ff", "darkorchid3",
           "magenta", "blue4","coral1","cornflowerblue",
           "#984ea3","#ff7f00","#ffff33","#f781bf","green",
           "deeppink","darkblue","cyan","deepskyblue", "darkgoldenrod")
############################################################################

# PCoA plot with symbol control in ggplot2
# plot by site_name
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=10))+
  geom_point(aes(color=site_name), size = 4) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme(plot.title = element_text(size = rel(2)),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=12),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig
# plot by rootstock
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=8))+
  geom_point(aes(color=rootstock), alpha = 0.7, size = 4) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig
# plot by cultivar
fig <- ggplot(data=map2, aes(x=ax1.scores, y=ax2.scores, size=8))+
  geom_point(aes(color=cultivar)) + 
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")+
  theme_bw(base_size=12)+
  ggtitle("16S rRNA PCoA (Bray-Curtis)")
fig

# PERMANOVA
adonis(otu_dist~map2$site_name)
adonis(otu_dist~map2$rootstock)
adonis(otu_dist~map2$cultivar)

# MAKE A BARPLOT USING PHYLOSEQ
# INSTALL PACKAGES
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)


# ADD THE TAXONOMY W/O PREFIX
# WRITE.CSV
write.csv(taxonomy, file = "16S_taxonomy.csv")
# READ .CSV
tax_16S = read.csv("16S_tax.csv", sep=',', header=T)
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
# IMPORT OTU TABLE
otu_phyl <- import_biom("/Users/arifinabintarti/Documents/Parent/apple_replant/16S analysis/single_rare20.biom")
otu_phyl
# ADD MAP
rownames(map20) <- map20$sample_code
phyloseq_map <- sample_data(map20)
otu_map <- merge_phyloseq(OTU,TAX,phyloseq_map)
otu_map
otu_map.ra <- transform_sample_counts(otu_map, function(x) x/sum(x))

# CHANGE THE RANK NAMES
rank_names(otu_map)
colnames(tax_table(otu_map)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
sample_phylum <- tax_glom(otu_map, taxrank = "Phylum",NArm = FALSE)
# MERGE SAMPLES BY SITE_NAMES
site_phyl <- merge_samples(otu_map, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(otu_map)$site_name)
# MERGE TAXA BY RANK
site_phylum <- tax_glom(site_phyl, taxrank = "Phylum", NArm = FALSE)
site_RA <- transform_sample_counts(site_phylum, function(x) 100 * x/sum(x))
# MAKE A BAR PLOT
plot_bar(site_RA, "site_name", fill="Phylum")+
  theme(axis.text.y=element_text(size=22),axis.text.x=element_text(size=22,angle=49,hjust =1), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=16),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))

#### BOXPLOT #####
library(ggplot2)
library(plyr)
library(phyloseq)

otu_map

# get abundance in %
otu_RA <- transform_sample_counts(otu_map, function(x) x/sum(x))
site_phyl <- merge_samples(otu_RA, "site_name")
sample_data(site_phyl)$site_name <- levels(sample_data(otu_map)$site_name)
# agglomerate taxa
glom <- tax_glom(site_phyl, taxrank = 'Phylum', NArm = FALSE)

# create dataframe from phyloseq object
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)

# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~Phylum, function(x) c(median=median(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
remainder <- medians[medians$median <= 0.01,]$Phylum

# change their name to "Remainder"
dat[dat$Phylum %in% remainder,]$Phylum <- 'Remainder'

# boxplot
ggplot(dat,
       aes(x=Phylum,
           y=Abundance)) + geom_boxplot() + coord_flip()
















# CORE MICROBIOTA ANALYSIS 
# INSTALL PACKAGES
library(BiocInstaller)
source("http://www.bioconductor.org/biocLite.R")
useDevel()
biocLite("microbiome")
library(microbiome)

# OCCUPANCY VS ABUNDANCE
# comulative Occ_Abund
otu_PA <- 1*((otu>0)==1)
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Occ <- rowSums(otu_PA)/ncol(otu_PA)
otu <- otu[rowSums(otu)>0,]
otu_rel <- decostand(otu, method="total", MARGIN=2)
com_abund <- rowSums(otu_rel)

#color code for the top most abundant OTUs
color_top <- com_abund
color_top[] <- 'black' 
top_abund <- log10(com_abund)[order((com_abund), decreasing=TRUE)][1:10] 
color_top[names(color_top) %in% names(top_abund)] <- 'red'
plot(log10(com_abund), Occ, col=color_top, pch=20, ylab='Occupancy', xlab='log(sum relative abundace per OTU)\n', 
     main='occupancy-abundance\n plot (comulative)')

#Mean Occ_Abund
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Mean_otu_PA <- apply(otu_PA, 1, mean)
Mean_Occ <- rowSums(otu_PA)/ncol(otu_PA)
Mean_abund <- apply(otu_rel, 1, mean)



############################################################################
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
#otu abundances and frequencies required for generating the basic plots
otu_map.ra = transform_sample_counts(otu_map, function(x) x/sum(x))
otu_table(otu_map.ra)
otu_map.ra

#phylum.sum = tapply(taxa_sums(site_phyl), tax_table(site_phyl)[, "Phylum"], sum, NArm=FALSE)
#top15phyla = names(sort(phylum.sum, TRUE))[1:15]

abund_val <- function(normalized){
  otu.abun = apply(otu_table(normalized),1,mean)
  otu.freq = rowSums(otu_table(normalized) != 0)/45
  phyla = as.vector(data.frame(tax_table(normalized))$Phylum)
  levels(phyla) = c(levels(phyla),"Other")
  keephyla = c("Bacteroidetes","Proteobacteria","Actinobacteria", "Acidobacteria", "Thaumarchaeota", "Gemmatimonadetes", "Latescibacteria", "Planctomycetes", "Verrucomicrobia", "Chloroflexi", "Nitrospirae")
  phyla[!(phyla %in% keephyla)] = "Other"
  phyla = as.vector(phyla)
  phyla=as.factor(phyla)
  otuabun = cbind.data.frame(abundance=log(otu.abun),frequency=otu.freq,phyla)
  return(otuabun)
}

abund_all <- abund_val(otu_map.ra)
# Use color brewer to pick a color scheme for the phyla
brew = brewer.pal(12, "Paired")
# Create a scatterplot of OTUs showing their average relative abundance and occupancy
ggplot(abund_all, aes(x=abundance,y=frequency,color=phyla)) + 
  geom_point(size=3) + xlab("Mean relative abundance (log10 scale)") + 
  ylab("Mean occupancy (n=45)") + scale_colour_brewer(palette="Paired")+
  labs(title="Bacterial abundancy vs. occupancy plots")+ xlim(-15, -4.6)+
  theme(plot.title = element_text(size = rel(2)),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
      legend.text=element_text(size=18),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))











#########################################################################################################
# calculate general relative abundances of Phylum
otu_phyl= tax_glom(otu_map, "Phylum", NArm = FALSE)
tax_table(otu_phyl) 
otu_phyl_RA = transform_sample_counts(otu_phyl, function(x) 100 * x/sum(x)) 
RA=taxa_sums(otu_phyl_RA)/45 
sort(RA, decreasing = TRUE)
sample_sums(otu_phyl_RA)
tax_table(otu_phyl_RA)
# Relative Abundance per OTU (F17, F28, F31, F39, F32, F41, F07, F30, F29, F40) --> 10 most abundant OTU
dat <- psmelt(otu_map.ra_Phyl)
dat
dat$Phylum <- as.character(dat$Phylum)
Phylum_abundance <- aggregate(Abundance~Sample+Phylum, dat, FUN=sum)
Phylum_abundance
library(reshape)
Phylum_abundance <- cast(Phylum_abundance, Sample~Phylum)
###########################################################################################################

# PROCRUSTES - PROTEST
# Classical MDS 16S
otu_dist <- vegdist(t(otu), method='bray')
otu_pcoa <- cmdscale(otu_dist, eig=T)
otu_pcoa
# Classical MDS ITS
otu_distITS <- vegdist(t(otuITS), method='bray')
otu_pcoaITS <- cmdscale(otu_distITS, eig=T)
otu_pcoaITS
# Procrustest
pro <- procrustes(otu_pcoa, otu_pcoaITS, scale = TRUE, symmetric = TRUE)
summary(pro)
plot(pro)

protest=protest(otu_pcoa, otu_pcoaITS, scale=TRUE, symmetric=TRUE, permutations=999)
plot(protest)


# PCoA PLOT USING PHYLOSEQ AND GGPLOT
plot_ordination(
  physeq = otu_map,
  ordination = otu_pcoa,
  color = "site_name",
  title = "16S rRNA PCoA (Bray-Curtis)"
) +
  geom_point(aes(color = site_name), size = 5) +
  geom_point(colour = "grey90", size = 1.5) +
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  scale_size(guide="none")







# BIOGEOGRAPHICAL ANALYSIS
ordinate <- read.table('LONG_LAT.csv', sep=',', header=TRUE)
test_ordinate <- read.table('test.csv', sep=',', header=TRUE)
test_ordinate
ordinate1 <- ordinate[,1:4]
library(ggmap)
library(maps)

apporchard_map <- get_map(location = c(lon = -85.99,lat = 43.57),zoom = 14,
                          source="google", maptype = "satellite", crop=FALSE, scale=2)
ggmap(apporchard_map)
dev.off()
bbox <- make_bbox(ordinate1$LON, ordinate1$LAT, f = 0.24)
map <- get_map(bbox)
ggmap(map)
ggmap(map) +
  geom_point(data = ordinate1, aes(x = LON, y = LAT, color=site_name, alpha = 0.8), size = 8, shape = 20) +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=15),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))+
  guides(fill=FALSE, alpha=FALSE, size=FALSE)



# MANTEL TEST
library('ade4')
spatial.dists <- dist(cbind(ordinate1$LON, ordinate1$LAT))
spatial.dists
otu_dist
otu_distITS
mantel.rtest(spatial.dists, otu_dist, nrepet = 999)
mantel.rtest(spatial.dists, otu_distITS, nrepet = 999)


######## Occupancy-Abundancy #########

#Comulative Occ_Abund
otu_PA <- 1*((otu>0)==1)
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Occ_otu <- rowSums(otu_PA)/ncol(otu_PA)
abs_otu <- otu[rowSums(otu)>0,]
otu_rel <- decostand(abs_otu, method="total", MARGIN=2)
com_abund_otu <- rowSums(otu_rel)

#color code for the top most abundant OTUs
color_top <- com_abund_otu
color_top[] <- 'black' 
otu_top_abund <- log10(com_abund_otu)[order((com_abund_otu), decreasing=TRUE)][1:10] #top 10 most abundant OTUs 
color_top[names(color_top) %in% names(otu_top_abund)] <- 'red'
plot(log10(com_abund_otu), Occ_otu, col=color_top, pch=20, ylab='Occupancy (n=45)', xlab='log(sum relative abundace per OTU)\n (n=20838 OTUs)', 
     main='Bacterial and archaeal occupancy-abundance\n plot (comulative)')
#Mean Occ_Abund
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Mean_otu_PA <- apply(otu_PA, 1, mean)
Mean_Occ <- rowSums(otu_PA)/ncol(otu_PA)
Mean_abund <- apply(otu_rel, 1, mean)

#order(Mean_abund[names(Mean_abund) %in% unique(selected_otus_switch$otu)])

#Creating df for plotting with ggplot and adding color code for the shared and unique OTUs
df_occ <- data.frame(otu=names(Occ_otu), occ=Occ_otu) 
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
postscript('bacterial and archaeal_occ_abund.eps', width = 4.5, height = 4)
brew = brewer.pal(12, "Paired")
ggplot(data=occ_abun, aes(x=abun, y=occ, color=phyla)) +
  geom_point(size=4, pch=20)+
  scale_colour_brewer(palette="Paired")+
  labs(x="log(mean relative abundace per OTU)\n (n=20838 OTUs)", y= "Mean occupancy (n=45)", title='Bacterial abundancy vs. occupancy plots') +
  geom_text_repel(data=occ_abun[occ_abun$col=='red',], aes(label=otu), box.padding = unit(0.45, "lines"), show.legend = FALSE) +
  geom_hline(aes(yintercept=.4), linetype='dashed', size=1.5) +
  geom_vline(aes(xintercept=-2.5), linetype='dashed', size=1.5) + 
  theme(plot.title = element_text(size = rel(2)),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"))+
  scale_y_continuous(breaks=seq(0,1,.2)) +
  scale_x_continuous(breaks=seq(-10,2,2))
dev.off()



abund_val <- function(normalized){
  otu.abun = apply(otu_table(normalized),1,mean)
  otu.freq = rowSums(otu_table(normalized) != 0)/45
  phyla = as.vector(data.frame(tax_table(normalized))$Phylum)
  levels(phyla) = c(levels(phyla),"Other")
  keephyla = c("Bacteroidetes","Proteobacteria","Actinobacteria", "Acidobacteria", "Thaumarchaeota", "Gemmatimonadetes", "Latescibacteria", "Planctomycetes", "Verrucomicrobia", "Chloroflexi", "Nitrospirae")
  phyla[!(phyla %in% keephyla)] = "Other"
  phyla = as.vector(phyla)
  phyla=as.factor(phyla)
  otuabun = cbind.data.frame(abundance=log(otu.abun),frequency=otu.freq,phyla)
  return(otuabun)
}
abund_all <- abund_val(otu_map.ra)
# Use color brewer to pick a color scheme for the phyla
brew = brewer.pal(12, "Paired")
# Create a scatterplot of OTUs showing their average relative abundance and occupancy
ggplot(abund_all, aes(x=abundance,y=frequency,color=phyla)) + 
  geom_point(size=3) + xlab("Mean relative abundance (log10 scale)") + 
  ylab("Mean occupancy (n=45)") + scale_colour_brewer(palette="Paired")+
  labs(title="Bacterial abundancy vs. occupancy plots")+ xlim(-15, -4.6)+
  theme(plot.title = element_text(size = rel(2)),axis.text=element_text(size=22), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=18),legend.title = element_text(size = 14), legend.spacing.x = unit(1.0, 'cm'))

abund_all
phyla = as.vector(data.frame(tax_table(otu_map.ra))$Phylum)
phyla
levels(phyla) = c(levels(phyla),"Other")
keephyla = c("Bacteroidetes","Proteobacteria","Actinobacteria", "Acidobacteria", "Thaumarchaeota", "Gemmatimonadetes", "Latescibacteria", "Planctomycetes", "Verrucomicrobia", "Chloroflexi", "Nitrospirae")
phyla[!(phyla %in% keephyla)] = "Other"
phyla = as.vector(phyla)
phyla=as.factor(phyla)

occ_abun <- cbind.data.frame(occ_abun, phyla, by='otu')
occ_abun





















