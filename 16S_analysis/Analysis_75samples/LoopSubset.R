#Example of a loop for Fina
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
source("https://bioconductor.org/biocLite.R")
biocLite()
install.packages('reshape')
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)
setwd('/Users/arifinabintarti/Documents/Parent/apple_replant/')
wd <- print(getwd())
otu <- read.table('otu_table_v1.txt', sep='\t', header=T, row.names = 1)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
head(otu)
dim(otu)
taxonomy <- otu[,'taxonomy']
taxonomy
otu <- otu[,-46]
tax=cbind(row.names(otu), as.character(taxonomy))
dim(tax)
#finding OTU IDs
#Top10=c("GQ871301.1.964", "FQ659953.1.1390", "HE614852.1.1516", "GQ302575.1.1520", "FM209322.1.1496", "FN421989.1.1411", "LN614619.1.1387"
#        tax[grep("LN614619.1.1387", tax[,1]),]
      
        
        set.seed(13)
        
        #summary(otu)
        #min. number of sequences is ~26K
        #subsample to an even number of sequences per sample
        otu_even28k <- t(rrarefy(t(otu), min(colSums(otu))))
        head(otu_even28k)
        sort(colSums (otu_even28k, na.rm = FALSE, dims = 1), decreasing = FALSE)
        #summary(otu_even26k)
        
        otu_even28k <- otu_even28k[rowSums(otu_even28k)>0,]
        otu_rare_PA <- 1*(otu_even28k>0)
        s <- specnumber(otu_even28k, MARGIN = 2)
        h <- diversity(t(otu_even28k), index = 'shannon')
        pielou <- h/log(s)
        
        map.div <- map

        dim(map)
map2 <- map[31:75,]
map3 <- map[c(56:57, 66:67),]
map4 <- map[c(62:65, 68:70),]
block=unique(map2[,"site_name"])

for (i in 1:length(block)){
  #subset the OTU table to only include those blocks
  temp=otu_even28k[,map2[,"site_name"]==block[i]]

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

