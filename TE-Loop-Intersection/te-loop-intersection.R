## ChIA-pet/TE intersection

# Load up data tables
loop_data = read.table("/nfs/boylelab_turbo/mouseENCODE.new/CTCF-TE-Analysis/chia_pet.dat", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(loop_data) = c("id_cp", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "iab", "fdr", "factor_cp", "species", "cell")

loop_te_data = read.table("/nfs/boylelabnr_turbo/ChIA-PET/ctcf_loop_te_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(loop_te_data) = c("id_rp", "chrom", "chromstart_rp", "chromend_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromstart_ctcf", "chromend_ctcf", "signalvalue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr")

loop_motif_data = read.table("/nfs/boylelab_turbo/mouseENCODE.new/CTCF-TE-Analysis/ctcf_loop_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(loop_motif_data) = c("id_cp", "chrom", "start1", "end1", "start2", "end2", "iab", "fdr", "cell", "species", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "id_mt", "start_mt", "end_mt", "score_mt", "strand_mt")

loop_te_motif_data = read.table("/nfs/boylelab_turbo/mouseENCODE.new/CTCF-TE-Analysis/ctcf_loop_te_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(loop_te_motif_data) = c("id_rp", "chrom", "chromStart_rp", "chromEnd_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr", "id_mt", "chromStart_mt", "chromEnd_mt", "strand_mt", "score_mt")

# Create GM12878-loop-centric data table for TE-loop-motif intersections
source("build_extended_loop_datatable.R")
dat.GM12878 = add_te_intersections(loop_data[which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"),c(1:4,6:7)], loop_te_data)
dat.GM12878 = add_motif_orientations(dat.GM12878, loop_motif_data)

# Generate the pie chart for Figure 2.
pdf("TE-derived-loop-fractions.GM12878.pdf")
pie(c(nrow(dat.GM12878[which(dat.GM12878$te_left & dat.GM12878$te_right),]), nrow(dat.GM12878[which((dat.GM12878$te_left | dat.GM12878$te_right) & !(dat.GM12878$te_left & dat.GM12878$te_right)),]), nrow(dat.GM12878[which(!dat.GM12878$te_left & !dat.GM12878$te_right),])), labels = c("TE Both", "TE One", "No TE"), main = "Fraction of GM12878 Loops")
dev.off()

# Generate tabular data regarding motif convergence for Figure 2 A, B, and C (Also shown in supplemental table 6).
write.table(get_motif_data_table(dat.GM12878), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE, file="TE-derived-loop-fraction_motif-layouts.GM12878.txt")

# Repeat tabular data compilation for K562.
dat.K562 = add_te_intersections(loop_data[which(loop_data$cell == "K562" & loop_data$factor == "RAD21"),c(1:4,6:7)], loop_te_data)
dat.K562 = add_motif_orientations(dat.K562, loop_motif_data)
write.table(get_motif_data_table(dat.K562), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE, file="TE-derived-loop-fraction_motif-layouts.K562.txt")


## Hi-C loop/TE intersection for CH12

# Read in data
hic_te_data = read.table("/nfs/boylelab_turbo/mouseENCODE.new/CTCF-TE-Analysis/ctcf_hic_te_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(hic_te_data) = c("id_rp", "chrom", "chromstart_rp", "chromend_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromstart_ctcf", "chromend_ctcf", "signalvalue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr")

hic_motif_data = read.table("/nfs/boylelab_turbo/mouseENCODE.new/CTCF-TE-Analysis/ctcf_hic_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(hic_motif_data) = c("id_cp", "chrom", "start1", "end1", "start2", "end2", "iab", "fdr", "cell", "species", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "id_mt", "start_mt", "end_mt", "score_mt", "strand_mt")

hic_te_motif_data = read.table("/nfs/boylelab_turbo/mouseENCODE.new/CTCF-TE-Analysis/ctcf_hic_te_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(hic_te_motif_data) = c("id_rp", "chrom", "chromStart_rp", "chromEnd_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr", "id_mt", "chromStart_mt", "chromEnd_mt", "strand_mt", "score_mt")

# Compile CH12 tabular data for supplemental table 6.
dat.CH12 = add_te_intersections(loop_data[which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C"),c(1:4,6:7)], loop_te_data)
dat.CH12 = add_motif_orientations(dat.CH12, loop_motif_data)
write.table(get_motif_data_table(dat.CH12), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE, file="TE-derived-loop-fraction_motif-layouts.CH12.txt")

# Write data tables to files for use in loop conservation analysis.
write.table(dat.CH12, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, file="data.CH12.txt")
write.table(dat.GM12878, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, file="data.GM12878.txt")
write.table(dat.K562, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, file="data.K562.txt")
