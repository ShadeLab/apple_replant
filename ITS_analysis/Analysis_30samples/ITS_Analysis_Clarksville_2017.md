Raw Sequences Analysis of ITS gene from Apple Rhizosphere Soil in Clarksville
There are 30 samples (C01-C30). The raw sequences of the samples along with 45 samples (F01-F45) (75 soil samples in total) are stored in this directory:
/mnt/research/ShadeLab/Sequence/raw_sequence/AppleReplant/20180328_ITS_PE/

# Analysis of ITS Miseq Data

1. all the analyses results are stored in:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_30samples_Clarksville

2. primer/adapter removal was performed using cutadapt 1.17 with Python 2.7.15.
3. trimming was conducted with length of 200 bp
4. OTU table was built using both open reference of UNITE V.7.2 and de novo clustering with 97% of similarity
5. Classifying was conducted using "consensus_taxonomy.txt" file produced by CONSTAX
6. further analyses were conducted using QIIME (v.1.9.1) platform.
7. unfortunately, taxonomy file (consensus_taxonomy.txt) produced from CONSTAX is biom unfriendly so it cannot be added into OTU table biom file using command: biom add-metadata
8. add taxonomy to OTU table can be conducted/modified using R software, Excel, and Phyloseq 

# QUALITY CHECKING AND PRE-FILTERING

a) make a directory call 'rawreads' and put the unzipped fastq file in there
b) count read numbers
```
for fastq in rawreads/*.fastq
do wc -l $fastq
done > reads_raw.counts
```
c) produce reads quality graphs using FastQC
```
mkdir stats

cat rawreads/*R1_001.fastq > raw_reads_R1.fastq; cat rawreads/*R2_001.fastq > raw_reads_R2.fastq

module load FastQC/0.11.7-Java-1.8.0_162

fastqc raw_reads_R1.fastq raw_reads_R2.fastq -o stats && rm -rf raw_reads_R1.fastq raw_reads_R2.fastq
```

## MERGE PAIRED READS 
```
mkdir mergedfastq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs rawreads/*R1*.fastq -reverse rawreads/*R2*.fastq -fastq_maxdiffs 10 -fastq_minmergelen 50 -relabel @ -fastqout mergedfastq/merged.fastq
```
### output
```
4823385  Pairs (4.8M)
   3417605  Merged (3.4M, 70.85%)
   1840814  Alignments with zero diffs (38.16%)
   1255531  Too many diffs (> 10) (26.03%)
    123219  No alignment found (2.55%)
         0  Alignment too short (< 16) (0.00%)
     27030  Merged too short (< 50)
    417024  Staggered pairs (8.65%) merged & trimmed
    207.89  Mean alignment length
    283.49  Mean merged length
      0.96  Mean fwd expected errors
      0.74  Mean rev expected errors
      0.22  Mean merged expected errors
```
## CHECK SEQUENCE QUALITY OF MERGED SEQS USING USEARCH AND VSEARCH
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged.fastq -output stats_eestats2_USEARCH.txt

/mnt/home/bintarti/vsearch-2.8.4-linux-x86_64/bin/vsearch -fastq_stats mergedfastq/merged2.fastq -log stats_results_VSEARCH.txt
```
### output
```
4601240 reads, max len 484, avg 281.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4555392( 99.0%)    4597024( 99.9%)    4601216(100.0%)
   100    4272714( 92.9%)    4421917( 96.1%)    4440687( 96.5%)
   150    4117291( 89.5%)    4371651( 95.0%)    4438364( 96.5%)
   200    3944832( 85.7%)    4245254( 92.3%)    4360838( 94.8%)
   250    3770636( 81.9%)    4091828( 88.9%)    4238271( 92.1%)
   300     677838( 14.7%)     794362( 17.3%)     867429( 18.9%)
   350     173782(  3.8%)     218912(  4.8%)     256250(  5.6%)
   400      15972(  0.3%)      22710(  0.5%)      29544(  0.6%)
   450       1319(  0.0%)       2299(  0.0%)       3547(  0.1%)
```
## REMOVE PRIMER AND ADAPTERS WITH CUTADAPT, CHECK THE STATS
```
#upgrade pip
pip install --upgrade pip

#install latest version of cutadapt 2.0
pip3 install --user --upgrade cutadapt
```
```
CS1-ITS1 (fwd): 5’- CTTGGTCATTTAGAGGAAGTAA – 3’ (EMP/Smith and Peay 2014)
CS2-ITS2 (rev): 5’- GCTGCGTTCTTCATCGATGC – 3’ (EMP/Smith and Peay 2014)
Reverse complement of adapter-reverse primer: GCATCGATGAAGAACGCAGC

/mnt/home/bintarti/.local/bin/cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o cut_merged.fastq mergedfastq/merged.fastq > cut_adpt_results.txt

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 cut_merged.fastq -output cutdapt_eestats2_USEARCH.txt

#output: 4600167 reads, max len 463, avg 239.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4384372( 95.3%)    4435805( 96.4%)    4440797( 96.5%)
   100    4243581( 92.2%)    4411100( 95.9%)    4438871( 96.5%)
   150    4046648( 88.0%)    4289781( 93.3%)    4367229( 94.9%)
   200    3862516( 84.0%)    4138727( 90.0%)    4256775( 92.5%)
   250    1450277( 31.5%)    1606508( 34.9%)    1693344( 36.8%)
   300     203478(  4.4%)     248866(  5.4%)     286097(  6.2%)
   350      39477(  0.9%)      52476(  1.1%)      64854(  1.4%)
   400       3066(  0.1%)       4801(  0.1%)       6918(  0.2%)
   450          0(  0.0%)          1(  0.0%)          1(  0.0%)
```

## FILTERING, TRIMMING, QUALITY CHECK
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter cut_merged.fastq -fastq_maxee 1 -fastq_trunclen 200 -fastq_maxns 0 -fastaout filtered_cut_merged.fa -fastqout filtered_cut_merged.fastq

###output:
   100.0% Filtering, 89.2% passed
   3417190  Reads (3.4M)                    
    282509  Discarded reads length < 200
        34  Discarded read with > 0 Ns
     86167  Discarded reads with expected errs > 1.00
   3048480  Filtered reads (3.0M, 89.2%)

#quality check
fastqc filtered_cut_merged.fastq
```
# CLUSTERING AND DENOISING (OTUs)

## DEREPLICATE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_cut_merged.fastq -fastaout derep_filtered_cut_merged.fasta -sizeout

#output: 3048480 seqs, 431421 uniques, 290795 singletons (67.4%)
Min size 1, median 1, max 351199, avg 5.47
```
## REMOVE SINGLETONS
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize derep_filtered_cut_merged.fasta -fastaout nosingl_derep_filtered_cut_merged.fasta -minsize 2

#output: Sorting 140626 sequences
```
## OPEN REFERENCE-BASED OTU PICKING AGAINST FUNGAL ITS DATABASE UNITE (v.8.0) at 97% IDENTITY
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global derep_filtered_cut_merged.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v.8.0/sh_refs_qiime_ver8_97_02.02.2019.fasta  -strand plus -uc ref_seqs.uc -dbmatched UNITE_reference.fasta -notmatched UNITE_failed_closed.fq

# output: 100.0% Searching, 56.7% matched
```
### SORT BY SIZE AND CLUSTERING-MAKING OTU TABLE at 97% SEQUENCE SIMILARITY (de novo otu-picking)-UPARSE ON SEQUENCE THAT FAILED TO HIT REFERENCE/UNITE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -sortbysize UNITE_failed_closed.fq -fastaout sorted_UNITE_failed_closed.fq

# output: Sorting 186735 sequences

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus  sorted_UNITE_failed_closed.fq -minsize 2 -otus DENOVO_otus.fasta -uparseout uparse_otus.txt -relabel OTU_

#output: 100.0% 1890 OTUs, 327 chimeras
```
## Combine the rep sets between de novo and reference-based OTU picking
```
cat UNITE_reference.fasta DENOVO_otus.fasta > REPRESENTATIVE_SET.fna

/mnt/home/bintarti/python_scripts-master/fasta_number.py REPRESENTATIVE_SET.fna OTU_ > NUMBERED_REP_SET.fasta
```
## CONSTRUCT OTU TABLE (MAP READS BACK TO REPRESENTATIVE SEQUENCES AS DATABASE FILE)
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged.fastq -db NUMBERED_REP_SET.fasta -strand plus -id 0.97 -uc OTU_map.uc -otutabout OPEN_REF_OTU_TABLE_ITS.txt

#output: 92.8% matched
3170055 / 3417605 mapped to OTUs (92.8%) 
```
## TAXONOMIC CLASSIFICATION USING CONSTAX
```
please refer to how "Running CONSTAX on the MSU HPCC on lab guru : https://my.labguru.com/knowledge/documents/330

# add "NUMBERED_REP_SET.fasta" into: mnt/home/bintarti/CONSTAX_hpcc/otus
# change the PATH in: mnt/home/bintarti/CONSTAX_hpcc/config
# output directories: outputs, taxonomy_assignments, training_files
# move the output directories to: /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Workflow2_ITS1_OpenReferenceClust/CONSTAX_OTU
# Use the taxonomy file : outputs/consensus_taxonomy.txt
```
## CONVERT .TXT FILE INTO .BIOM FILE - CONSTAX output is not BIOM friendly, thus do not use this command: biom add_metadata
```
biom convert -i OPEN_REF_OTU_TABLE_ITS.txt -o OPEN_REF_OTU_TABLE_ITS.biom --table-type="OTU table" --to-json

#output file: OPEN_REF_OTU_TABLE_ITS.biom
```
## ALIGN SEQUENCES USING MUSCLE 

install muscle3.8.31_i86linux64 in your home directory
```
/mnt/home/bintarti/muscle3.8.31_i86linux64 -in NUMBERED_REP_SET.fasta -out OTUS_aligned.fasta -maxiters 2 -diags1
```
# START THE ANALYSES USING QIIME 1.9.1 PIPELINE

First, you should install QIIME 1.9.1 using Python2.7 environment

##activate qiime1 environment 
```
conda activate qiime1
```
## FILTER EXCESS GAPS FROM ALIGNMENT USING QIIME 1.9.1
```
filter_alignment.py -i OTUS_aligned.fasta -o OTUS_filtered_alignment
```
## MAKE PHYLOGENY WITH FASTREE
```
make_phylogeny.py -i OTUS_filtered_alignment/OTUS_aligned_pfiltered.fasta -o OTUS_rep_set.tre
```
## Rarefy OTU table to lowest sequencing depth
```
# summarize OTU table
biom summarize-table -i OPEN_REF_OTU_TABLE_ITS.biom -o OTU_table_sum.txt

# rarefaction
single_rarefaction.py -d 59874 -o OTU_ITS_rarefied.biom -i OPEN_REF_OTU_TABLE_ITS.biom

# summarize rarefied OTU table
biom summarize-table -i OTU_ITS_rarefied.biom -o OTU_ITS_rarefied_sum.txt

# Convert rarefied OTU table file: .biom into .txt for further analyses on R
biom convert -i OTU_ITS_rarefied.biom -o OTU_ITS_rarefied.txt --to-tsv

# Convert OPEN_REF_OTU_TABLE_ITS.biom file into .txt for further analyses on R
biom convert -i OPEN_REF_OTU_TABLE_ITS.biom -o rarefied/OPEN_REF_OTU_TABLE_ITS.txt --to-tsv
```