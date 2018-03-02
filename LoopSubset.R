#Example of a loop for Fina
install.packages(c('vegan', 'tidyverse'))
<<<<<<< HEAD
install.packages('reshape')
source("https://bioconductor.org/biocLite.R")
biocLite()
=======
#install.packages('reshape')
#source("https://bioconductor.org/biocLite.R")
#biocLite()
>>>>>>> 328757c190f948356fffac8b754a48bf793674be
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

<<<<<<< HEAD
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/')
=======
#setwd('/Users/arifinabintarti/Documents/apple_replant/')
>>>>>>> 328757c190f948356fffac8b754a48bf793674be
otu <- read.table('otu_table_v2.txt', sep='\t', header=T, row.names = 1)
map <- read.table('clean_map_data.csv', header=TRUE, sep=",")
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
tax=cbind(row.names(otu), as.character(taxonomy))
#finding OTU IDs
#Top10=c("GQ871301.1.964", "FQ659953.1.1390", "HE614852.1.1516", "GQ302575.1.1520", "FM209322.1.1496", "FN421989.1.1411", "LN614619.1.1387"
#        tax[grep("LN614619.1.1387", tax[,1]),]
        
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
<<<<<<< HEAD
        dim(map)
map2 <- map[c(1:55, 58:65, 68:75),] 
map3 <- map[c(56:57, 66:67),]
block=unique(map3[,"site_name"])

for (i in 1:length(block)){
  #subset the OTU table to only include those blocks
  temp=otu_even26k[,map3[,"site_name"]==block[i]]
=======

block=unique(map[,"site_name"])

for (i in 1:length(block)){
  #subset the OTU table to only include those blocks
  temp=otu_even26k[,map[,"site_name"]==block[i]]
>>>>>>> 328757c190f948356fffac8b754a48bf793674be
  
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

