## ChIA-pet/TE intersection, as presented in Fig. 4, Sup. Fig. 3, and Sup. Table 5.
require(ggplot2)
require(reshape)

# Load Functions
source("build_extended_loop_datatable.R")

# Load up data tables
loop_data = read.table("../data/hadoop/chia_pet.dat", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_cp", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "iab", "fdr", "factor_cp", "species", "cell"))
loop_te_data = read.table("../data/hadoop/results/ctcf_loop_te_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_rp", "chrom", "chromstart_rp", "chromend_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromstart_ctcf", "chromend_ctcf", "signalvalue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr"))
loop_motif_data = read.table("../data/hadoop/results/ctcf_loop_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_cp", "chrom", "start1", "end1", "start2", "end2", "iab", "fdr", "cell", "species", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "id_mt", "start_mt", "end_mt", "score_mt", "strand_mt"))
loop_te_motif_data = read.table("../data/hadoop/results/ctcf_loop_te_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_rp", "chrom", "chromStart_rp", "chromEnd_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr", "id_mt", "chromStart_mt", "chromEnd_mt", "strand_mt", "score_mt"))
hic_te_data = read.table("../data/hadoop/results/ctcf_hic_te_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_rp", "chrom", "chromstart_rp", "chromend_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromstart_ctcf", "chromend_ctcf", "signalvalue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr"))
hic_motif_data = read.table("../data/hadoop/results/ctcf_hic_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_cp", "chrom", "start1", "end1", "start2", "end2", "iab", "fdr", "cell", "species", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "id_mt", "start_mt", "end_mt", "score_mt", "strand_mt"))
hic_te_motif_data = read.table("../data/hadoop/results/ctcf_hic_te_motif_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_rp", "chrom", "chromStart_rp", "chromEnd_rp", "name_rp", "class_rp", "family_rp", "pctdiv_rp", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "factor", "species", "cell", "id_cp", "start1", "end1", "start2", "end2", "iab", "fdr", "id_mt", "chromStart_mt", "chromEnd_mt", "strand_mt", "score_mt"))
enr_te_fams = read.table("enriched_te_fams.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, col.names=c("species", "name"), row.names="name")

# Create cell-wise data tables
dat.GM12878 = add_te_intersections(loop_data[which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"),c(1:4,6:7)], loop_te_data)
dat.GM12878 = add_motif_orientations(dat.GM12878, loop_motif_data)
dat.K562 = add_te_intersections(loop_data[which(loop_data$cell == "K562" & loop_data$factor == "RAD21"),c(1:4,6:7)], loop_te_data)
dat.K562 = add_motif_orientations(dat.K562, loop_motif_data)
dat.CH12 = add_te_intersections(loop_data[which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C"),c(1:4,6:7)], loop_te_data)
dat.CH12 = add_motif_orientations(dat.CH12, loop_motif_data)

# Create single data table with uniform column conventions
tmp1 = add_ctcf_summits(dat.CH12, loop_ctcf_data, col="id_cp")
tmp1$species = "mm9"
tmp1$cell_q = "CH12"
tmp1$cell_t1 = "GM12878"
tmp1$cell_t2 = "K562"
tmp1 = tmp1[,c(1:10, 34:37, 11:12, 21:22 ,13:20, 23:33)]
cn = colnames(tmp1)
cn[15:26] = c("l_maps_t1", "r_maps_t1", "l_maps_t2", "r_maps_t2", "l_overlaps_t1", "r_overlaps_t1", "l_overlaps_t2", "r_overlaps_t2", "class_t1", "class_t2", "row_t1", "row_t2")
colnames(tmp1) = cn

tmp2 = add_ctcf_summits(dat.K562, loop_ctcf_data, col="id_cp")
tmp2$species = "hg19"
tmp2$cell_q = "K562"
tmp2$cell_t1 = "CH12"
tmp2$cell_t2 = "GM12878"
tmp2$r_maps_t2 = tmp2$l_maps_t2 = TRUE
tmp2 = tmp2[,c(1:10, 32:35, 11:12, 36:37, 13:20, 21:31)]
cn = colnames(tmp2)
cn[15:26] = c("l_maps_t1", "r_maps_t1", "l_maps_t2", "r_maps_t2", "l_overlaps_t1", "r_overlaps_t1", "l_overlaps_t2", "r_overlaps_t2", "class_t1", "class_t2", "row_t1", "row_t2")
colnames(tmp2) = cn

tmp3 = add_ctcf_summits(dat.GM12878, loop_ctcf_data, col="id_cp")
tmp3$species = "hg19"
tmp3$cell_q = "GM12878"
tmp3$cell_t1 = "CH12"
tmp3$cell_t2 = "K562"
tmp3$r_maps_t2 = tmp3$l_maps_t2 = TRUE
tmp3 = tmp3[,c(1:10, 32:35, 13:14, 36:37, 15:18, 11:12, 19:20, 21:31)]
cn = colnames(tmp3)
cn[15:26] = c("l_maps_t1", "r_maps_t1", "l_maps_t2", "r_maps_t2", "l_overlaps_t1", "r_overlaps_t1", "l_overlaps_t2", "r_overlaps_t2", "class_t1", "class_t2", "row_t1", "row_t2")
colnames(tmp3) = cn

annotated_loop_data = rbind(tmp1, tmp2, tmp3)

# Draw pie chart for Fig. 4A-C.
pdf("TE-derived-loop-fractions.pdf")
pie(c(nrow(annotated_loop_data[which(annotated_loop_data$te_left & annotated_loop_data$te_right),]), nrow(annotated_loop_data[which((annotated_loop_data$te_left | annotated_loop_data$te_right) & !(annotated_loop_data$te_left & annotated_loop_data$te_right)),]), nrow(annotated_loop_data[which(!annotated_loop_data$te_left & !annotated_loop_data$te_right),])), labels = c("TE Both", "TE One", "No TE"), main = "Fraction of Loops")
dev.off()
# Convergence data and loop counts for motif convergence in Figs. 4A-C.
tmp = get_motif_data_table(annotated_loop_data)
tmp$pct_cnv = tmp$Convergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_tnd = tmp$Tandem / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_div = tmp$Divergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pc = tmp[,4] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pd = tmp[,5] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_unk = tmp[,6] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
write.table(tmp, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE, file="motif-orientations.all.txt")

# Draw pie chart for Supplementary Fig. 3A.
pdf("TE-derived-loop-fractions.GM12878.pdf")
pie(c(nrow(dat.GM12878[which(dat.GM12878$te_left & dat.GM12878$te_right),]), nrow(dat.GM12878[which((dat.GM12878$te_left | dat.GM12878$te_right) & !(dat.GM12878$te_left & dat.GM12878$te_right)),]), nrow(dat.GM12878[which(!dat.GM12878$te_left & !dat.GM12878$te_right),])), labels = c("TE Both", "TE One", "No TE"), main = "Fraction of GM12878 Loops")
dev.off()
# Loop counts and convergence data for Supplementary Fig. 3A and Sup Table 5.
tmp = get_motif_data_table(dat.GM12878)
tmp$pct_cnv = tmp$Convergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_tnd = tmp$Tandem / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_div = tmp$Divergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pc = tmp[,4] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pd = tmp[,5] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_unk = tmp[,6] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
write.table(tmp, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE, file="motif-orientations.GM12878.txt")

# Draw pie chart for Supplementary Fig. 3B.
pdf("TE-derived-loop-fractions.K562.pdf")
pie(c(nrow(dat.K562[which(dat.K562$te_left & dat.K562$te_right),]), nrow(dat.K562[which((dat.K562$te_left | dat.K562$te_right) & !(dat.K562$te_left & dat.K562$te_right)),]), nrow(dat.K562[which(!dat.K562$te_left & !dat.K562$te_right),])), labels = c("TE Both", "TE One", "No TE"), main = "Fraction of Loops")
dev.off()
# Loop counts and convergence data for Supplementary Fig. 3B and Sup Table 5.
tmp = get_motif_data_table(dat.K562)
tmp$pct_cnv = tmp$Convergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_tnd = tmp$Tandem / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_div = tmp$Divergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pc = tmp[,4] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pd = tmp[,5] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_unk = tmp[,6] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
write.table(tmp, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE, file="motif-orientations.K562.txt")

# Draw pie chart for Supplementary Fig. 3C.
pdf("TE-derived-loop-fractions.CH12.pdf")
pie(c(nrow(dat.CH12[which(dat.CH12$te_left & dat.CH12$te_right),]), nrow(dat.CH12[which((dat.CH12$te_left | dat.CH12$te_right) & !(dat.CH12$te_left & dat.CH12$te_right)),]), nrow(dat.CH12[which(!dat.CH12$te_left & !dat.CH12$te_right),])), labels = c("TE Both", "TE One", "No TE"), main = "Fraction of Loops")
dev.off()
# Loop counts and convergence data for Supplementary Fig. 3C and Sup Table 5.
tmp = get_motif_data_table(dat.CH12)
tmp$pct_cnv = tmp$Convergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_tnd = tmp$Tandem / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_div = tmp$Divergent / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pc = tmp[,4] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_pd = tmp[,5] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
tmp$pct_unk = tmp[,6] / (tmp$Convergent + tmp$Tandem + tmp$Divergent + tmp[,4] + tmp[,5] + tmp[,6])
write.table(tmp, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE, file="motif-orientations.CH12.txt")

# Stacked bars for Fig. 4D.
tmp.CH12 = table(c(dat.CH12$te_name_l, dat.CH12$te_name_r))
tmp.CH12 = tmp.CH12[order(tmp.CH12)]
tmp.CH12 = tmp.CH12[which(names(tmp.CH12) %in% enr_te_fams$name)]
tmp.GM12878 = table(c(dat.GM12878$te_name_l, dat.GM12878$te_name_r))
tmp.GM12878 = tmp.GM12878[order(tmp.GM12878)]
tmp.GM12878 = tmp.GM12878[which(names(tmp.GM12878) %in% enr_te_fams$name)]
tmp.K562 = table(c(dat.K562$te_name_l, dat.K562$te_name_r))
tmp.K562 = tmp.K562[order(tmp.K562)]
tmp.K562 = tmp.K562[which(names(tmp.K562) %in% enr_te_fams$name)]
dat = enr_te_fams
dat$K562 = dat$GM12878 = dat$CH12 = NA
colnames(dat)[3:5] = c("CH12","GM12878","K562")
rownames(dat) = unique(c(names(tmp.CH12), names(tmp.GM12878), names(tmp.K562)))
dat[names(tmp.CH12),"CH12"] = tmp.CH12
dat[names(tmp.GM12878),"GM12878"] = tmp.GM12878
dat[names(tmp.K562),"K562"] = tmp.K562
tmp = melt(dat)
tmp$name = factor(tmp$name, levels=dat$name)

pdf("loop_te-enriched-species-intersections.pdf")
ggplot(tmp, aes(x=variable, y=value, fill=species)) +
geom_bar(stat="identity", position="fill")
dev.off()

# Save objects used elsewhere
save(dat.GM12878, dat.CH12, dat.K562, annotated_loop_data, file="loop_intersection.Rdata")