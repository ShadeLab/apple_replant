# Analysis of fungal ITS-Apple replant samples (KearnsPat)

All analyses were performed in FASTX toolkit (v. 0.0.14), ITSx (v. 1.0.11), and USEARCH v. 10.0.240 (64-bit) with the UPARSE OTU picking method. Subsequent analyses were performed in QIIME (v. 1.8).

raw sequence data stored on HPCC

apple replant soil ITS:  45 total files (F01-F30)


Copied and renamed 45 apple replant ITS sequences to working space for analysis
ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Pats_workflow/

## Merge Paired End Reads
```
#decompress the reads
gunzip *.gz

mkdir mergedfastq

#make the usearch program executable
ln -s /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 usearch10

#will use a more variable sequence size 200-500bp due to the variability of ITS amplicons.
./usearch10 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 200 -fastq_maxmergelen 500
```
### output
```
6818052  Pairs (6.8M)
   4373305  Merged (4.4M, 64.14%)
   2563437  Alignments with zero diffs (37.60%)
   1854285  Too many diffs (> 10) (27.20%)
    234054  No alignment found (3.43%)
         0  Alignment too short (< 16) (0.00%)
    356408  Merged too short (< 200)
         0  Merged too long (> 500)
         0  Exp.errs. too high (max=1.0) (0.00%)
    578828  Staggered pairs (8.49%) merged & trimmed
    208.28  Mean alignment length
    290.58  Mean merged length
      0.90  Mean fwd expected errors
      0.69  Mean rev expected errors
      0.24  Mean merged expected errors
```
## Quality filter fastq sequences and convert to fasta format with the FASTX toolkit (v. 0.0.14)
```
# module load FASTX/0.0.14
fastq_quality_filter -q 30 -i mergedfastq/merged.fq -o mergedfastq/merged_fil.fq -p 50
```
## Dereplicate sequences
```
./usearch10 -fastx_uniques mergedfastq/merged_fil.fasta -fastaout mergedfastq/derep_merged_fil.fasta -sizeout

#output: 4373162 seqs, 1031920 uniques, 774826 singletons (75.1%)
```
## Remove Singeltons
```
./usearch10 -sortbysize mergedfastq/derep_merged_fil.fasta -fastaout mergedfastq/nosigs_derep_merged_fil.fasta -minsize 2

#output: Sorting 257094 sequences
```
## Precluster Sequences
```
./usearch10 -cluster_fast mergedfastq/nosigs_derep_merged_fil.fasta -centroids mergedfastq/denoised_nosigs_derep_merged_fil.fasta -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

#output: 257094 seqs (tot.size 3598336), 257094 uniques, 0 singletons (0.0%)
Seqs  257094 (257.1k)
  Clusters  10240 (10.2k)
  Max size  609712 (609.7k)
  Avg size  351.4
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  157Mb
      Time  30.0s
Throughput  8569.8 seqs/sec.
```
## Reference-based OTU picking against UNITE fungal ITS database (v. 7.2) at 97% sequence similarity
```
./usearch10 -usearch_global mergedfastq/denoised_nosigs_derep_merged_fil.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v7.2/sh_refs_qiime_ver7_97_s_01.12.2017.fasta -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatched mergedfastq/failed_closed.fq
```
## Sort by size and then de novo OTU picking on sequences that failed to hit UNITE
```
./usearch10 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

./usearch10 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

#output: 2710 OTUs, 512 chimeras
```
## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```
## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
./usearch10 -usearch_global mergedfastq/merged.fq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

#count OTUs
grep -c '>' mergedfastq/full_rep_set.fna
#OTUs: 3511
```
# Switch to QIIME (v. 1.9.1)

## Assign taxonomy to UNITE (v. 7.2) and UCLUST
```
assign_taxonomy.py -i mergedfastq/full_rep_set.fna -o taxonomy -r /mnt/research/ShadeLab/UNITE_v7.2/sh_refs_qiime_ver7_97_s_01.12.2017.fasta -t /mnt/research/ShadeLab/UNITE_v7.2/sh_taxonomy_qiime_ver7_97_s_01.12.2017.txt

#count Unassigned taxa
grep -c 'Unassigned' taxonomy/full_rep_set_tax_assignments.txt
#2692
```
## Add taxonomy to OTU table
```
biom convert -i OTU_table.txt -o otu_table_from_txt_json.biom --table-type="OTU table" --to-json

biom add-metadata -i otu_table_from_txt_json.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```
## Filter Rhizaria
```
filter_taxa_from_otu_table.py -i otu_table_tax.biom -o otu_table_tax_filt2.biom -n k__Rhizaria,Unassigned

# summarize otu table:
biom summarize-table -i otu_table_tax_filt2.biom -o otu_table_sum2.txt
```
## Rarefy OTU table to lowest sequencing depth
```
single_rarefaction.py -d 56728 -o single_rare.biom -i otu_table_tax_filt.biom
```
## Convert single_rare.biom into single_rare_its.txt
```
biom convert -i single_rare.biom -o OTU_single_rare.txt --to-tsv --header-key taxonomy
```
