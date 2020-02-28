#!/bin/bash

# CTCF motif-word enrichment pipeline.

DATA_DIR=../../data

# First, make CTCF motif predictions within CTCF-bound transposon insertions.

bedtools intersect -a <(zcat $DATA_DIR/repeatMasker/hg19.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk '{split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\n", $1, $2, $3, A[1]}}') -b <(awk '{printf "%s\t%d\t%d\n", $1, ($2+$10)-1, ($2+$10)+1}' $DATA_DIR/ChIP-seq/CTCF/hg19.merged.narrowPeak) -wa -u > ctcf_bound_repeats.hg19.bed

bedtools intersect -a <(zcat $DATA_DIR/repeatMasker/mm9.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk '{split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\n", $1, $2, $3, A[1]}}') -b <(awk '{printf "%s\t%d\t%d\n", $1, ($2+$10)-1, ($2+$10)+1}' $DATA_DIR/ChIP-seq/CTCF/mm9.merged.narrowPeak) -wa -u > ctcf_bound_repeats.mm9.bed

bedtools getfasta -fi $DATA_DIR/motifs/FIMO/hg19.fa -bed ctcf_bound_repeats.hg19.bed -name -fo hg19.bound_rmsk.fa
bedtools getfasta -fi $DATA_DIR/motifs/FIMO/mm9.fa -bed ctcf_bound_repeats.mm9.bed -name -fo mm9.bound_rmsk.fa

mkdir {bound_rmsk_hg19,bound_rmsk_mm9}
WD=$(pwd)

cd bound_rmsk_hg19
/usr/bin/env fimo --max-stored-scores 1000000 $DATA_DIR/motifs/CTCF.ren.meme $WD/hg19.bound_rmsk.fa
cd $WD

cd bound_rmsk_mm9
/usr/bin/env fimo --max-stored-scores 1000000 $DATA_DIR/motifs/CTCF.ren.meme $WD/mm9.bound_rmsk.fa
cd $WD

cat bound_rmsk_hg19/fimo_out/fimo.txt | awk '{if (NR > 1) {printf "%s\t%d\t%d\t.\t%s\t%s\t%s\n", $2, $3, $4, $7, $5, $9}}' > CTCF.bound_rmsk.hg19.bed
cat bound_rmsk_mm9/fimo_out/fimo.txt | awk '{if (NR > 1) {printf "%s\t%d\t%d\t.\t%s\t%s\t%s\n", $2, $3, $4, $7, $5, $9}}' > CTCF.bound_rmsk.mm9.bed

# Next, normalize counts, calculate odds ratios, and run
# Fisher's exact tests for enrichment in TE families within R.
/usr/bin/env Rscript --vanilla calc-word-enrichments.R

# That's it!
