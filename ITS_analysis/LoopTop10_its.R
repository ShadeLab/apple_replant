# LOOP

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
tax=cbind(row.names(otu), as.character(taxonomy))
dim(tax)
#finding OTU IDs
#Top10=c("GQ871301.1.964", "FQ659953.1.1390", "HE614852.1.1516", "GQ302575.1.1520", "FM209322.1.1496", "FN421989.1.1411", "LN614619.1.1387"
#        tax[grep("LN614619.1.1387", tax[,1]),]


set.seed(13)

#summary(otu)
#subsample to an even number of sequences per sample
head(sort(rowSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(colSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
otu_rare_PA <- 1*(otu>0)
s <- specnumber(otu, MARGIN = 2)
h <- diversity(t(otu), index = 'shannon')
pielou <- h/log(s)

map.div <- map

map
map1 <- map[1:43,]
map2 <- map[c(31:38),]

map3 <- map[c(50:53),]
map4 <- map[c(46:47),]
map5 <- map[c(54:55),]
map6 <- map[c(56:57),]
map7 <- map[c(68:75),]
block=unique(map[,"site_name"])

for (i in 1:length(block)){
  #subset the OTU table to only include those blocks
  temp=otu[,map[,"site_name"]==block[i]]
  
  #how many samples in this block?
  samples=ncol(temp)
  
  #how many total sequences in the whole dataset?
  totalseq=sum(temp)
  
  #what is the number of sequences per OTU?
  r=rowSums(temp)
  
  #sort the otus by abundance, with the most abundant first
  s=sort(r, decreasing=TRUE, index.return=TRUE)
  
  #using the returned index from the sort (s$ix), sort the dataset with the most abundant OTU first
  temp.s=temp[s$ix,]
  
  #pull out the top ten
  temp.10=temp.s[1:10,]
  
  #what is the cumulative contribution of the top ten OTUs?
  cumulcontri=sum(temp.10)/totalseq
  print(cumulcontri)
  
  #who are these top 10 OTUs?
  names=row.names(temp.10)
  id=NULL
  for(j in 1:length(names)){
    id.temp=tax[grep(names[j],tax[,1]),]
    id=rbind(id, id.temp)
  }
  colnames(id)=c("otuID", "SilvaTaxonomy")
  write.table(id, file= paste(block[i],".txt", sep=""),sep="\t", quote=FALSE)
  
  #vector for output data
  out=c(samples,cumulcontri)
  
}

