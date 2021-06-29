####
# Pipeline to assemble orthology data for CTCF sites used to
# generate Figure 2A.

####
# First step is finding the intersections between CTCF sites
# in each species, and with TEs...

# Copy in data and add unique identifiers to each row.
DATA_HG=../data/ChIP-seq/CTCF/hg19.merged.narrowPeak
DATA_MM=../data/ChIP-seq/CTCF/mm9.merged.narrowPeak

awk 'BEGIN {N=1}; {printf "%s\t%d\t%d\t%d\t%d\t.\t%f\t%f\t%f\t%d\n", $1, $2, $3, N, $5, $7, $8, $9, $10, $11; N++}' $DATA_HG > hg19.merged.narrowPeak
awk -v N=$(tail -n 1 dat.hg | awk '{print $4}') 'BEGIN {N++}; {printf "%s\t%d\t%d\t%d\t%d\t.\t%f\t%f\t%f\t%d\n", $1, $2, $3, N, $5, $7, $8, $9, $10, $11; N++}' $DATA_MM > mm9.merged.narrowPeak

# Cross-map CTCF peaks with bnMapper. The version used has been modified from the
# original, and is available at https://github.com/Boyle-Lab/bx-python.
# liftOver chains were prepared as described in ../data/liftover/README.
# Mouse data were mapped to human using the set of swapped human to mouse
# chains to come as close as possible to reversible mapping between species.
bnMapper.py CTCF.hg19.merged.narrowPeak ../data/liftover/hg19.mm9.rbest.chain.gz -f narrowPeak -k -i narrowPeak > CTCF.hg19.merged.lifted.bed
bnMapper.py CTCF.mm9.merged.narrowPeak ../data/liftover/hg19.mm9.rbest.rev.chain.gz -f narrowPeak -k -i narrowPeak > CTCF.mm9.merged.lifted.bed

# Merge mapped peaks with native peaks to get union sets
bedtools merge -i <(cat CTCF.hg19.merged.narrowPeak CTCF.mm9.merged.lifted.bed | awk '{printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, $2+$10}' | bedtools sort) -c 4,5 -o collapse -delim ","  | awk '{psum = 0; N=split($5, A, ","); for (i=0; i<=N; i++) {psum += A[i]}; printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, (psum/N)-$2}' >  CTCF.hg19.mm9.union.bed
bedtools merge -i <(cat CTCF.mm9.merged.narrowPeak CTCF.hg19.merged.lifted.bed | awk '{printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, $2+$10}' | bedtools sort) -c 4,5 -o collapse -delim "," | awk '{psum = 0; N=split($5, A, ","); for (i=0; i<=N; i++) {psum += A[i]}; printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, (psum/N)-$2}'  >  CTCF.mm9.hg19.union.bed

# Label columns to indicate species-specific/shared occupancy
# Human-referenced
bedtools intersect -v -a CTCF.hg19.merged.narrowPeak -b CTCF.mm9.merged.lifted.bed | awk '{printf "%s\t%d\t%d\t%d\t%d\t1\t0\n", $1, $2, $3, $4, $10}' > tmp
bedtools intersect -v -b CTCF.hg19.merged.narrowPeak -a CTCF.mm9.merged.lifted.bed | awk '{printf "%s\t%d\t%d\t%d\t%d\t0\t1\n", $1, $2, $3, $4, $10}' >> tmp
bedtools intersect -v -a CTCF.hg19.mm9.union.bed -b tmp | awk '{printf "%s\t%d\t%d\t%d\t%d\t1\t1\n", $1, $2, $3, $4, $5}' > tmp2
cat tmp tmp2 > CTCF.hg19.mm9.union.labelled.bed
# Mouse-referenced
bedtools intersect -v -a CTCF.mm9.merged.narrowPeak -b CTCF.hg19.merged.lifted.bed | awk '{printf "%s\t%d\t%d\t%d\t%d\t1\t0\n", $1, $2, $3, $4, $10}' > tmp
bedtools intersect -v -b CTCF.mm9.merged.narrowPeak -a CTCF.hg19.merged.lifted.bed | awk '{printf "%s\t%d\t%d\t%d\t%d\t0\t1\n", $1, $2, $3, $4, $10}' >> tmp
bedtools intersect -v -a CTCF.mm9.hg19.union.bed -b tmp | awk '{printf "%s\t%d\t%d\t%d\t%d\t1\t1\n", $1, $2, $3, $4, $5}' > tmp2
cat tmp tmp2 > CTCF.mm9.hg19.union.labelled.bed

# Intersect with transposable elements
# Human-referenced
bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", $1, ($2+$5), ($2+$5)+1, $4, $5, $6, $7}' CTCF.hg19.mm9.union.labelled.bed) -b <(zcat ../data/repeatMasker/hg19.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wa -u > tmp
bedtools intersect -v -a CTCF.hg19.mm9.union.labelled.bed -b tmp | awk '{printf "%s\t0\n", $0}' > tmp2
bedtools intersect -a CTCF.hg19.mm9.union.labelled.bed -b tmp -wa -u | awk '{printf "%s\t1\n", $0}' >> tmp2
mv tmp2 CTCF.hg19.mm9.union.labelled.te.bed
# Mouse-referenced
bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", $1, ($2+$5), ($2+$5)+1, $4, $5, $6, $7}' CTCF.mm9.hg19.union.labelled.bed) -b <(zcat ..//data/repeatMasker/mm9.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wa -u > tmp
bedtools intersect -v -a CTCF.mm9.hg19.union.labelled.bed -b tmp | awk '{printf "%s\t0\n", $0}' > tmp2
bedtools intersect -a CTCF.mm9.hg19.union.labelled.bed -b tmp -wa -u | awk '{printf "%s\t1\n", $0}' >> tmp2
mv tmp2 CTCF.mm9.hg19.union.labelled.te.bed

# Add back unmapped features.
# Human-referenced
# Step 1: prepare a minimal data frame
awk '{printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, $10}' CTCF.mm9.merged.narrowPeak > tmp1
awk '{printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, $10}' CTCF.mm9.merged.lifted.bed > tmp2
# Step 2 is in R. See get_species-specfic_sites.R  -- output is written to tmp3
/usr/bin/env Rscript --vanilla get_species-specific_sites.R

# Step 3: Intersect with transposable elements.
bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", $1, ($2+$5), ($2+$5)+1, $4, $5, $6, $7}' tmp3) -b <(zcat ../data/repeatMasker/mm9.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wa -u > tmp
# Step 4: Label rows in tmp (mouse binding data) with TE overlap labels (tmp and tmp3 already have species occupancy labels)
bedtools intersect -v -a tmp3 -b tmp | awk '{printf "%s\t0\n", $0}' > tmp2
bedtools intersect -a tmp3 -b tmp -wa -u | awk '{printf "%s\t1\n", $0}' >> tmp2
# Step 5: Combine all outputs into complete data table and remove intermediates.
cat CTCF.hg19.mm9.union.labelled.te.bed tmp2 > tmp
mv tmp CTCF.hg19.mm9.union.labelled.te.bed
rm tmp*

# Mouse-referenced
# Step 1
awk '{printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, $10}' CTCF.hg19.merged.narrowPeak > tmp1
awk '{printf "%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $4, $10}' CTCF.hg19.merged.lifted.bed > tmp2
# Step 2, as above: Run commands in get_species-specfic_sites.R
/usr/bin/env Rscript --vanilla get_species-specific_sites.R

# Step 3
bedtools intersect -a <(awk '{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", $1, ($2+$5), ($2+$5), $4, $5, $6, $7}' tmp3) -b <(zcat ../data/repeatMasker/hg19.rmsk.bed.gz | grep -v "random" | grep -v "hap" | grep -v "Un" | grep -v "chrM" | awk 'BEGIN {i=1}; {split($4, A, "."); if (A[2] != "Low_complexity" && A[2] != "Satellite" && A[2] != "Simple_repeat" && A[2] != "tRNA" && A[2] != "rRNA" && A[2] != "scRNA" && A[2] != "snRNA" && A[2] != "srpRNA") {printf "%s\t%d\t%d\t%s\t%s\t%s\t%f\thg19\n", $1, $2, $3, A[1], A[2], A[3], $5; i++}}') -wa -u > tmp
# Step 4
bedtools intersect -v -a tmp3 -b tmp | awk '{printf "%s\t0\n", $0}' > tmp2
bedtools intersect -a tmp3 -b tmp -wa -u | awk '{printf "%s\t1\n", $0}' >> tmp2
# Step 5
cat CTCF.mm9.hg19.union.labelled.te.bed tmp2 > tmp
mv tmp CTCF.mm9.hg19.union.labelled.te.bed
rm tmp*

# Plot the Venn diagram -- commands are in plot_orthology-te_data.R.
# Manual adjustments to plot areas to make areas proportional to subset size were
# performed by hand in Adobe Illustrator to produce the final plot presented in Figre 1A.
/usr/bin/env Rscript --vanilla plot_orthology-te_data.R
