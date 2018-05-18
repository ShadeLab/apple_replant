# Analysis of ITS_apple_replant_adopted from KearnsPat

All analyses were performed in FASTX toolkit (v. 0.0.14), ITSx (v. 1.0.11), and USEARCH v. 10.0.240 (64-bit) with the UPARSE OTU picking method. Subsequent analyses were performed in QIIME (v. 1.8).

raw sequence data stored on HPCC

apple replant soil ITS:  75 total files


Copied apple replant ITS sequences to working space for analysis
ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_rawreads/


## Merge Paired End Reads
Here we are using the RDPs usearch64, and copied an executable to this directory.  All commands start with "./usearch64" instead of "usearch"

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

```
#decompress the reads
gunzip *.gzfastq_quality_filter -q 30 -i mergedfastq/merged.fq -o mergedfastq/merged_fil.fq -p 50

mkdir mergedfastq

#will use a more variable sequence size 200-500bp due to the varibility of ITS amplicons.
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs rawreads/*R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 200 -fastq_maxmergelen 500
```
### output

```
11641437  Pairs (11.6M)
   7720872  Merged (7.7M, 66.32%)
   4404251  Alignments with zero diffs (37.83%)
   3109816  Too many diffs (> 10) (26.71%)
    357273  No alignment found (3.07%)
         0  Alignment too short (< 16) (0.00%)
    453476  Merged too short (< 200)
         0  Merged too long (> 500)
         0  Exp.errs. too high (max=1.0) (0.00%)
    995852  Staggered pairs (8.55%) merged & trimmed
    208.90  Mean alignment length
    288.98  Mean merged length
      0.91  Mean fwd expected errors
      0.69  Mean rev expected errors
      0.23  Mean merged expected errors
```

## Quality filter fastq sequences and convert to fasta format with the FASTX toolkit (v. 0.0.14)

```
#module load FASTX/0.0.14
fastq_quality_filter -q 30 -i mergedfastq/merged.fq -o mergedfastq/merged_fil.fq -p 50

#convert fastq to fasta for ITSx (skip this step)
fastq_to_fasta -i mergedfqfastq/merged_fil.fq -o mergedfastq/merged_fa.fasta
```

## Use ITSx to remove non-ITS ribosomal fragments from fasta reads (skip this step)

```
#module load GNU/4.4.5
#module load ITSx/1.0.11
mkdir mergedfastq/ITSx
#this takes a while, so submit as job script
ITSx -i mergedfastq/merged_fa.fasta -o mergedfastq/ITSx/merged_fa_ITSx --cpu 16
```

## Dereplicate Sequence/find unique read sequences

```
 /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques mergedfastq/merged_fa.fasta -fastaout mergedfastq/uniques_combined_merged.fasta -sizeout
```
### output

```
00:33 3.1Gb   100.0% DF
00:33 3.2Gb  7720675 seqs, 1723759 uniques, 1282913 singletons (74.4%)
00:33 3.2Gb  Min size 1, median 1, max 633718, avg 4.48
00:38 3.2Gb   100.0% Writing mergedfastq/uniques_combined_merged.fasta
```

## Remove Singeltons

```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize mergedfastq/uniques_combined_merged.fasta -fastaout mergedfastq/nosigs_uniques_combined_merged.fasta -minsize 2
```

### output
```
00:05 669Mb   100.0% Reading mergedfastq/uniques_combined_merged.fasta
00:05 635Mb  Getting sizes                                            
00:05 649Mb  Sorting 440846 sequences
00:06 651Mb   100.0% Writing output
```
## Precluster Sequences

```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fasta -centroids mergedfastq/denoised_nosigs_uniques_combined_merged.fasta -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```
### output
```
00:03 194Mb  440846 seqs (tot.size 6437762), 440846 uniques, 0 singletons (0.0%)
Seqs  440846 (440.8k)
Clusters  14585 (14.6k)
Max size  1316697 (1.3M)
Avg size  441.4
Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
Max mem  255Mb Time  01:01
Throughput  7227.0 seqs/sec.
```
## Reference-based OTU picking against UNITE fungal ITS database (v. 7.2) at 97% sequence similarity

```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v7.2/sh_refs_qiime_ver7_97_s_01.12.2017.fasta -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatched mergedfastq/failed_closed.fq
```
## Sort by size and then de novo OTU picking on sequences that failed to hit GreenGenes

```
./mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

./mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
```
## Combine the rep sets between de novo and reference-based OTU picking

```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables

```
./mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged.fq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

#count OTUs
grep -c '>' mergedfastq/full_rep_set.fna
#OTUs
```

# Switch to QIIME (v. 1.8.0)

## Assign taxonomy to UNITE (v. 7.2) and UCLUST

```
assign_taxonomy.py -i mergedfastq/full_rep_set.fna -o taxonomy -r /mnt/research/ShadeLab/UNITE_v7.2/sh_refs_qiime_ver7_97_s_01.12.2017.fasta -t /mnt/research/ShadeLab/UNITE_v7.2/sh_taxonomy_qiime_ver7_97_s_01.12.2017.txt --uclust_similarity=0.5
#count Unassigned taxa
grep -c 'Unassigned' taxonomy/full_rep_set_tax_assignments.txt
#
```

## Add taxonomy to OTU table

```
biom convert -i OTU_table.txt --table-type='OTU table' -o otu_jsn.biom
biom add-metadata -i otu_jsn.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```

## Filter Protozoa and Unassigned taxa

```
filter_taxa_from_otu_table.py -i otu_table_tax.biom -o otu_table_tax_filt.biom -n k__Rhizaria,Unassigned

filter_taxa_from_otu_table.py -i otu_table_tax_filt.biom -o otu_table_tax_filt2.biom -n k__Protozoa,Unassigned

#summarize table
biom summarize-table -i otu_table_tax_filt2.biom -o otu_table_sum.txt
```

## Rarefy OTU table to lowest sequencing depth

```
single_rarefaction.py -d 5465 -o single_rare.biom -i otu_table_tax_filt2.biom
``` 

## Convert single_rare.biom into single_rare_its.txt

```
biom convert -i single_rare.biom -o single_rare_its.txt --header-key taxonomy -b
```

