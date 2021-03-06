# Analysis of 16S Miseq Data
raw sequence data stored on HPCC
/mnt/research/ShadeLab/Sequence/raw_sequence/AppleReplant/20171113_16S-V4_PE/

apple replant soil 16S DNA:  75 total files


Moved apple replant sequences to working space for analysis
ShadeLab/WorkingSpace/Bintarti/16S_analysis_apple_replant

All analysis results are in this path : 
/mnt/research/ShadeLab/WorkingSpace/Bintarti/16S_analysis_apple_replant



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
Totals:
4641655  Pairs (4.6M)
   3763326  Merged (3.8M, 81.08%)
   1380951  Alignments with zero diffs (29.75%)
    851888  Too many diffs (> 10) (18.35%)
     15877  No alignment found (0.34%)
         0  Alignment too short (< 16) (0.00%)
      9382  Merged too short (< 250)
      1182  Merged too long (> 300)
         0  Exp.errs. too high (max=1.0) (0.00%)
     11556  Staggered pairs (0.25%) merged & trimmed
    246.91  Mean alignment length
    253.09  Mean merged length
      0.50  Mean fwd expected errors
      1.38  Mean rev expected errors
      0.16  Mean merged expected errors
[bintarti@dev-intel14-phi my_apple_replant]$

```

## Dereplicate sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques mergedfastq/merged.fq -fastqout mergedfastq/uniques_combined_merged.fastq -sizeout
```
###Output
00:09 2.3Gb   100.0% Reading mergedfastq/merged.fq
00:16 2.5Gb   100.0% DF
00:16 2.5Gb  3763326 seqs, 1486603 uniques, 1094230 singletons (73.6%)
000:16 2.5Gb  Min size 1, median 1, max 12244, avg 2.53
00:41 2.5Gb   100.0% Writing mergedfastq/uniques_combined_merged.fastq


## Remove Singeltons
```
./usearch64 -sortbysize mergedfastq/uniques_combined_merged.fastq -fastqout mergedfastq/nosigs_uniques_combined_merged.fastq -minsize 2
```
### Output
00:04 938Mb   100.0% Reading mergedfastq/uniques_combined_merged.fastq
00:04 904Mb  Getting sizes
00:04 916Mb  Sorting 392373 sequences
00:04 918Mb   100.0% Writing output

## Precluster Sequences
```
./usearch64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```
### Output
Seqs  392373 (392.4k)
  Clusters  103889 (103.9k)
  Max size  23638 (23.6k)
  Avg size  25.7
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  505Mb
      Time  03:14
Throughput  2022.5 seqs/sec.


## Reference-based OTU picking using Silva database

```
./usearch64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db /mnt/research/ShadeLab/WorkingSpace/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta  -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatchedfq mergedfastq/failed_closed.fq
```
### Output
00:17 286Mb   100.0% Reading /mnt/research/ShadeLab/WorkingSpace/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fastaet97.fna
00:17 252Mb   100.0% Masking (fastnucleo)
00:27 253Mb   100.0% Word stats          
00:27 253Mb   100.0% Alloc rows
00:47 1.1Gb   100.0% Build index
00:47 1.1Gb  CPU has 28 cores, defaulting to 10 threads
06:47 1.2Gb   100.0% Searching denoised_nosigs_uniques_combined_merged.fastq, 37.1% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the Shade Lab working space.
Produce some output files - ref_seqs.uc (pre-clustered), closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking

## De novo OTU picking

```
#sort by size
./usearch64 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq
#cluster de novo
./usearch64 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
```
### Output
01:37 125Mb   100.0% 18778 OTUs, 27538 chimeras

## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
./usearch64 -usearch_global mergedfastq/merged.fq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom
```
outputs here = OTU table and OTU_jsn.biom -> for goin in qiime
00:01 66Mb    100.0% Reading mergedfastq/full_rep_set.fna
00:01 32Mb    100.0% Masking (fastnucleo)                
00:02 34Mb    100.0% Word stats          
00:02 34Mb    100.0% Alloc rows
00:03 108Mb   100.0% Build index
00:03 141Mb  CPU has 28 cores, defaulting to 10 threads
16:28 161Mb   100.0% Searching merged.fq, 92.7% matched
3487950 / 3763326 mapped to OTUs (92.7%)               
16:28 161Mb  Writing OTU_table.txt
16:28 161Mb  Writing OTU_table.txt ...done.
16:28 161Mb  Writing OTU_jsn.biom
16:28 161Mb  Writing OTU_jsn.biom ...done.


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
## Convert 
```
biom convert -i otu_table_tax_filt.biom -o otu_table_v2.txt --header-key taxonomy -b
```
## Align sequences to SILVA_128_QIIME_release with PyNAST 
```
align_seqs.py -i /mnt/research/ShadeLab/WorkingSpace/Bintarti/16S_analysis_apple_replant/mergedfastq/full_rep_set.fna -o alignment -t /mnt/home/bintarti/SILVA_128_QIIME_release/rep_set_aligned/97/97_otus_aligned.fasta
```

## Filter excess gaps from alignment
```
filter_alignment.py -i alignment/full_rep_set_aligned.fasta -o alignment/filtered_alignment
```

## Make phylogeny with fasttree
```
make_phylogeny.py -i alignment/filtered_alignment/full_rep_set_aligned_pfiltered.fasta -o rep_set.tre
```

## Summarize the OTU table and rarefy OTU table to lowest sequencing depth
```
biom summarize-table -i otu_table_tax_filt.biom -o otu_table_summary.txt

single_rarefaction.py -d 26235 -o single_rare.biom -i otu_table_tax_filt.biom
```

## Calculate alpha and beta diversity
```
beta_diversity.py -m bray_curtis,unweighted_unifrac,weighted_unifrac -i single_rare.biom -o beta_div -t rep_set.tre

alpha_diversity.py -m PD_whole_tree,shannon -i single_rare.biom -o alpha -t rep_set.tre
```

## Summarize taxonomy data
```
summarize_taxa.py -i otu_table_tax_filt.biom -o taxa_sum
```
## Creating a heatmap
```
make_otu_heatmap_html.py -i otu_table_tax_filt.biom -o OTU_Heatmap/
```
## Creating a taxa_chart
```
plot_taxa_summary.py -i taxa_sum/otu_table_tax_filt_L3.txt -l Phylum -o Taxa_Charts -k white
``
