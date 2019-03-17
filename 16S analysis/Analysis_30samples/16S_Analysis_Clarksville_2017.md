Raw Sequences Analysis of 16S rRNA gene from Apple Rhizosphere Soil in Clarksville
There are 30 samples (C01-C30). The raw sequences of the samples along with 45 samples (F01-F45) (75 soil samples in total) are stored in this directory:
/mnt/research/ShadeLab/Sequence/raw_sequence/AppleReplant/20171113_16S-V4_PE/

# Analysis of 16S Miseq Data

Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME1

Moved/copy apple replant sequences (30 soil samples (C01-C30) to working space for analysis:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/16S_30samples_Clarksville_2017

# Part I: Clustering

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

## 1) Merge Paired End Reads
```
# decompress the reads
gunzip *.gz

# make directory called "mergedfastq"
mkdir mergedfastq

# merge paired end reads
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -tabbedout mergedfastq/merged.report.txt -alnout mergedfastq/merged_aln.txt

# -fastq_maxdiffs 10: Allow 10 max differences in overlap region

### Output ###

 2227272  Pairs (2.2M)
   1822827  Merged (1.8M, 81.84%)
    686439  Alignments with zero diffs (30.82%)
    396405  Too many diffs (> 10) (17.80%)
      8040  No alignment found (0.36%)
         0  Alignment too short (< 16) (0.00%)
         0  Exp.errs. too high (max=1.0) (0.00%)
      4984  Staggered pairs (0.22%) merged & trimmed
    246.84  Mean alignment length
    253.03  Mean merged length
      0.51  Mean fwd expected errors
      1.34  Mean rev expected errors
      0.16  Mean merged expected errors
```
## 2) Check Sequence Quality of Merged Seqs
```
mkdir fastq_info
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged.fq -output fastq_info/eestats.txt

### output ###
1822827 reads, max len 473, avg 253.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    1798956( 98.7%)    1821060( 99.9%)    1822806(100.0%)
   100    1757791( 96.4%)    1811442( 99.4%)    1822464(100.0%)
   150    1727674( 94.8%)    1798124( 98.6%)    1820617( 99.9%)
   200    1703654( 93.5%)    1789076( 98.1%)    1819112( 99.8%)
   250    1657953( 91.0%)    1764501( 96.8%)    1811259( 99.4%)
   300        336(  0.0%)        433(  0.0%)        492(  0.0%)
   350         94(  0.0%)        144(  0.0%)        196(  0.0%)
   400          8(  0.0%)         32(  0.0%)         69(  0.0%)
   450          0(  0.0%)          0(  0.0%)          3(  0.0%)
```
## 3) Filter and Truncate the Merged Seqs
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter merged.fq -fastq_maxee 1 -fastq_trunclen 250 -fastqout filtered_merged.fq

### output ###
100.0% Filtering, 96.8% passed
   1822827  Reads (1.8M)                    
      4065  Discarded reads length < 250
     54261  Discarded reads with expected errs > 1.00
   1764501  Filtered reads (1.8M, 96.8%)
```
## 4) Dereplicate Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_merged.fq -fastqout uniques_filtered_merged.fastq -sizeout

### output ###
1764501 seqs, 678686 uniques, 481060 singletons (70.9%)
00:20 1.9Gb  Min size 1, median 1, max 6466, avg 2.60
```
## 5) Remove Singeltons
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize uniques_filtered_merged.fastq -fastqout nosigs_uniques_filtered_merged.fastq -minsize 2

### output ###
Sorting 197626 sequences
```
## 6) Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast nosigs_uniques_filtered_merged.fastq -centroids_fastq denoised_nosigs_uniques_filtered_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

### output ###
Seqs  197626 (197.6k)
  Clusters  61839 (61.8k)
  Max size  12059 (12.1k)
  Avg size  20.8
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  895Mb
      Time  48.0s
Throughput  4117.2 seqs/sec.
```
## 7) Closed Reference-based OTU Picking Using SILVA_132 Database
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global denoised_nosigs_uniques_filtered_merged.fastq -id 0.97 -db /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -strand plus -uc ref_seqs.uc -dbmatched SILVA_closed_reference.fasta -notmatchedfq failed_closed.fq

### output ###
43.4% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the home directory.
Produce some output files - ref_seqs.uc (pre-clustered), SILVA_closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking
```
## 8) De novo OTU picking
```
# sort by size
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize failed_closed.fq -fastaout sorted_failed_closed.fq

# cluster de novo
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus sorted_failed_closed.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

### output ###
11950 OTUs, 13960 chimeras
```
## 9) Combine the Rep Sets Between De novo and SILVA Reference-based OTU Picking
```
cat SILVA_closed_reference.fasta denovo_otus.fasta > FULL_REP_SET.fna
```
## 10) Map 'FULL_REP_SET.fna' Back to Pre-dereplicated Sequences and Make OTU Tables
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -usearch_global merged.fq -db FULL_REP_SET.fna -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

### output ###
1678493 / 1822827 mapped to OTUs (92.1%) 
```
# Part II: Switch to QIIME (v.1.9.1)

## 1) Assign taxonomy to SILVA_132_QIIME_release with PyNAST
```
assign_taxonomy.py -i FULL_REP_SET.fna -o taxonomy -r /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /mnt/home/bintarti/SILVA_132_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt
#-r reference -> path to silva db here (SILVA_128_QIIME_release)
# -t taxonomy
```
## 2) Convert OTU_table.txt. to OTU_table.from_txt_json.biom
```
biom convert -i OTU_table.txt -o OTU_table.biom --table-type="OTU table" --to-json
```
## 3) Add taxonomy to OTU table
```
biom add-metadata -i OTU_table.biom -o OTU_table_tax.biom --observation-metadata-fp=taxonomy/FULL_REP_SET_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```
## 4) Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i OTU_table_tax.biom -o OTU_table_tax_filt.biom -n D_4__Mitochondria,D_3__Chloroplast

# summarize OTU table

biom summarize-table -i OTU_table_tax.biom -o OTU_table_tax_sum.txt

biom summarize-table -i OTU_table_tax_filt.biom -o OTU_table_tax_filt_sum.txt
```
## 5) Rarefaction into the minimum number of sequence
```
single_rarefaction.py -d 26306 -o OTU_rarefied.biom -i OTU_table_tax_filt.biom

# summarize OTU table
biom summarize-table -i OTU_rarefied.biom -o OTU_rarefied_sum.txt
```
## 6) Convert and add taxonomy
```
biom convert -i OTU_rarefied.biom -o OTU_rarefied.txt --header-key taxonomy --to-tsv

biom convert -i OTU_table_tax_filt.biom -o OTU_table_tax_filt.txt --header-key taxonomy --to-tsv
```
## 7) Align sequences to SILVA_132_QIIME_release with PyNAST 
```
align_seqs.py -i FULL_REP_SET.fna -o alignment -t /mnt/home/bintarti/SILVA_132_QIIME_release/core_alignment/80_core_alignment.fna
```
## 8) Filter excess gaps from alignment
```
filter_alignment.py -i alignment/FULL_REP_SET_aligned.fasta -o alignment/filtered_alignment
```
## 9) Make phylogeny with fasttree
```
make_phylogeny.py -i alignment/filtered_alignment/FULL_REP_SET_aligned_pfiltered.fasta -o rep_set.tre
```
## 10) Summarize taxonomy data
```
summarize_taxa.py -i OTU_rarefied.biom -o taxa_sum
```





