#!/bin/bash

DATA_DIR=../../data

# Get all repeat names from the repeatMasker human and mouse annotations
zcat $DATA_DIR/repeatMasker/hg19.rmsk.bed.gz $DATA_DIR/repeatMasker/mm9.rmsk.bed.gz | awk '{split($4, A, "."); print A[1]}' | sort | uniq > reps.txt

zcat ../../data/repeatMasker/hg19.rmsk.bed.gz ../../data/repeatMasker/mm9.rmsk.bed.gz | awk '{split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") print A[1]}' | sort | uniq > reps.txt


# Extract ancestral sequences
extractAncestral.py <(zcat $DATA_DIR/repeatMasker/RMRBSeqs.embl) -s $(sort reps.txt | uniq | tr '\n' ',' | sed 's/,$//') > all_repeats_consensus.fa

# Score motifs in ancestral sequences.
score_motifs.pl $DATA_DIR/motifs/CTCF.ren.meme all_repeats_consensus.fa -scores -pseudo 0.001 -nthreads 8 -prefix "all-repeats_CTCF-ren"
# These will be further analyzed in R.
