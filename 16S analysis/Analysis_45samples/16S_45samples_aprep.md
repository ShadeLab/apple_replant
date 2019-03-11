16S_45samples_aprep

# Analysis of 16S Miseq Data
raw sequence data stored on HPCC
/mnt/research/ShadeLab/Sequence/raw_sequence/AppleReplant/20171113_16S-V4_PE/


Moved/copy apple replant sequences (45 soil samples (F01-F45)) to working space for analysis
ShadeLab/WorkingSpace/Bintarti/16S_45samples_aprep/20171113_16S-V4_PE/

# Part I: Clustering

## Merge Paired End Reads
Here we are using the RDPs usearch64, and copied an executable to this directory.  All commands start with "./usearch64" instead of "usearch"

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch


```
#decompress the reads
gunzip *.gz

mkdir mergedfastq
#rename files to remove -
rename "-" "" *fastq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64


./usearch64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 250 -fastq_maxmergelen 300

# call executable usearch64 from within this folder.
# merge pairs with R1 as the name.  Relabel with everything before the _ (samples will be named K11, K12 etc.) Allow 10 max differences in overlap region. put output file into the mergedfasq directory, named merged.fq. merged length between 250 & 300
```

### Output

```
 2414383  Pairs (2.4M)
   1945073  Merged (1.9M, 80.56%)
    694512  Alignments with zero diffs (28.77%)
    455483  Too many diffs (> 10) (18.87%)
      7837  No alignment found (0.32%)
         0  Alignment too short (< 16) (0.00%)
      5317  Merged too short (< 250)
       673  Merged too long (> 300)
         0  Exp.errs. too high (max=1.0) (0.00%)
      6572  Staggered pairs (0.27%) merged & trimmed
    246.90  Mean alignment length
    253.10  Mean merged length
      0.49  Mean fwd expected errors
      1.42  Mean rev expected errors
      0.16  Mean merged expected errors
```

## Dereplicate sequences

```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques mergedfastq/merged.fq -fastqout mergedfastq/uniques_combined_merged.fastq -sizeout
```
###Output

```
00:08 1.3Gb   100.0% DF
00:08 1.3Gb  1945073 seqs, 826701 uniques, 611974 singletons (74.0%)
00:08 1.3Gb  Min size 1, median 1, max 6734, avg 2.35
00:20 1.3Gb   100.0% Writing mergedfastq/uniques_combined_merged.fastq
```

## Remove Singeltons
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize mergedfastq/uniques_combined_merged.fastq -fastqout mergedfastq/nosigs_uniques_combined_merged.fastq -minsize 2
```
### Output

```
00:02 543Mb   100.0% Reading mergedfastq/uniques_combined_merged.fastq
00:02 509Mb  Getting sizes                                            
00:02 516Mb  Sorting 214727 sequences
00:03 517Mb   100.0% Writing output
```

## Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```
### Output
Seqs  214727 (214.7k)
  Clusters  69941 (69.9k)
  Max size  13056 (13.1k)
  Avg size  19.1
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  308Mb
      Time  01:06
Throughput  3253.4 seqs/sec.

## Reference-based OTU picking using Silva database

```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db /mnt/research/ShadeLab/WorkingSpace/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta  -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatchedfq mergedfastq/failed_closed.fq
```
### Output
00:05 284Mb   100.0% Reading /mnt/research/ShadeLab/WorkingSpace/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta
00:08 250Mb   100.0% Masking (fastnucleo)             
00:15 251Mb   100.0% Word stats          
00:15 251Mb   100.0% Alloc rows
00:27 1.1Gb   100.0% Build index
03:58 1.2Gb   100.0% Searching denoised_nosigs_uniques_combined_merged.fastq, 45.0% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the Shade Lab working space.
Produce some output files - ref_seqs.uc (pre-clustered), closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking

## De novo OTU picking

```
#sort by size
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

#cluster de novo
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
```
### Output
00:49 95Mb    100.0% 13858 OTUs, 12980 chimeras

## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -usearch_global mergedfastq/merged.fq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom
```
```
00:00 61Mb    100.0% Reading mergedfastq/full_rep_set.fna
00:01 27Mb    100.0% Masking (fastnucleo)                
00:01 28Mb    100.0% Word stats          
00:01 28Mb    100.0% Alloc rows
00:02 91Mb    100.0% Build index
00:02 124Mb  CPU has 28 cores, defaulting to 10 threads
08:34 137Mb   100.0% Searching merged.fq, 92.8% matched
1805245 / 1945073 mapped to OTUs (92.8%)               
08:34 137Mb  Writing OTU_table.txt
08:34 137Mb  Writing OTU_table.txt ...done.
08:34 137Mb  Writing OTU_jsn.biom
08:34 137Mb  Writing OTU_jsn.biom ...done.
```

# Part II: Switch to QIIME

## Assign taxonomy to SILVA_128_QIIME_release
```
assign_taxonomy.py -i full_rep_set.fna -o taxonomy -r /mnt/home/bintarti/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t /mnt/home/bintarti/SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt
#-r reference -> path to silva db here (SILVA_128_QIIME_release)
# -t taxonomy
```

## Add taxonomy to OTU table
## Convert OTU_table.txt. to OTU_table.from_txt_json.biom

```
biom convert -i OTU_table.txt -o OTU_table.from_txt_json.biom --table-type="OTU table"

biom add-metadata -i OTU_table.from_txt_json.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```
## Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i otu_table_tax.biom -o otu_table_tax_filt.biom -n o__Streptophyta,o__Chlorophyta,f__mitochondria,Unassigned
```
## Rarefaction
```
single_rarefaction.py -d 28172 -o single_rare20.biom -i otu_table_tax_filt.biom
```
## Convert and add taxonomy
```
biom convert -i single_rare20.biom -o otu_table_20.txt --header-key taxonomy -b
```

## Align sequences to SILVA_128_QIIME_release with PyNAST 
```
align_seqs.py -i /mnt/research/ShadeLab/WorkingSpace/Bintarti/16S_analysis_apple_replant/mergedfastq/full_rep_set.fna -o alignment -t /mnt/home/bintarti/SILVA_128_QIIME_release/rep_set_aligned/97/97_otus_aligned.fasta
```

## Filter excess gaps from alignment
```
filter_alignment.py -i alignment/full_rep_set_aligned.fasta -o alignment/filtered_alignment
```







