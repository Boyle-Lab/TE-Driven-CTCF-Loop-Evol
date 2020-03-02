#!/bin/bash

# Predict phylogenetic gain and loss of repeat elements in the human
# and mouse genomes.

DATA_DIR=../../data

# First pull out all genomic repeats overlapping CTCF sites in human
# and mouse genomes.
bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%f\t%d\n", $1, ($2+$10), ($2+$10)+1, $4, $5, $6, $7, $8, $9, $10}' $DATA_DIR/ChIP-seq/CTCF/hg19.merged.narrowPeak) -b <(zcat $DATA_DIR/repeatMasker/hg19.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wo | awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%f\t%d\n", $1, $2-25, $2+25, $4, $5, $6, $7, $8, $9, 25}' > tmp.hg19

bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%f\t%d\n", $1, ($2+$10), ($2+$10)+1, $4, $5, $6, $7, $8, $9, $10}' $DATA_DIR/ChIP-seq/CTCF/mm9.merged.narrowPeak) -b <(zcat $DATA_DIR/repeatMasker/mm9.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wo | awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%f\t%d\n", $1, $2-25, $2+25, $4, $5, $6, $7, $8, $9, 25}' > tmp.mm9

# Label TE-derived CTCF sites as orthologs, gains, and losses.
/usr/bin/envh mapGL.py tmp.hg19 "(((hg19:1,mm9:3):2,(canFam2:4,equCab2:5):6):90,loxAfr3:7)" hg19 mm9 $DATA_DIR/liftover/hg19.mm9.rbest.chain.gz $DATA_DIR/liftover/hg19.equCab2.rbest.chain.gz $DATA_DIR/liftover/hg19.canFam2.rbest.chain.gz $DATA_DIR/liftover/hg19.loxAfr3.rbest.chain.gz -i narrowPeak > hg19.gain-loss.out

/usr/bin/env mapGL.py tmp.mm9 "(((hg19:1,mm9:3):2,(canFam2:4,equCab2:5):6):90,loxAfr3:7)" mm9 hg19 $DATA_DIR/liftover/mm9.hg19.rbest-rev-hg19.chain.gz $DATA_DIR/liftover/mm9.equCab2.rbest.chain.gz $DATA_DIR/liftover/mm9.canFam2.rbest.chain.gz $DATA_DIR/liftover/mm9.loxAfr3.rbest.chain.gz -i narrowPeak > mm9.gain-loss.out

# Delete intermediate files
rm tmp.hg19 tmp.mm9

# Assemble output into a table for R
bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%d\t%d\t%s\n", $1, ($2+$5), ($2+$5)+1, $4, $5, $6}' hg19.gain-loss.out) -b <(zcat $DATA_DIR/repeatMasker/hg19.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wo | awk '{printf "%s\t%d\t%d\t%d\t%d\t%s\t1\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, $4, $5, $6, $10, $11, $12, $13}' > gain-loss_TE.txt

bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%d\t%d\t%s\n", $1, ($2+$5), ($2+$5)+1, $4, $5, $6}' mm9.gain-loss.out) -b <(zcat $DATA_DIR/repeatMasker/mm9.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\tmm9\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wo | awk '{printf "%s\t%d\t%d\t%d\t%d\t%s\t1\t%s\t%s\t%s\t%f\tmm9\n", $1, $2, $3, $4, $5, $6, $10, $11, $12, $13}' >> gain-loss_TE.txt

