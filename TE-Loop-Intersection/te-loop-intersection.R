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


## Add intersection with enriched TE types (Fig 2D, sup fig. 3)

# Read in data
enr_te_fams = read.table("enriched_te_fams.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
colnames(enr_te_fams) = c("species", "name")
rownames(enr_te_fams) = enr_te_fams$name

# Load commands
source("build_extended_loop_datatable.R")

# CH12
tmp = dat.CH12

# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"GM12878_l_maps"] == FALSE) | (tmp$te_right & tmp[,"GM12878_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp
tmp = add_te_fam(tmp, hic_te_data)
tmp = add_te_spec(tmp, enr_te_fams)
dat.CH12 = tmp

# GM12878
tmp = dat.GM12878
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp
tmp = add_te_fam(tmp, loop_te_data)
tmp = add_te_spec(tmp, enr_te_fams)
dat.GM12878 = tmp

# K562
tmp = dat.K562
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp
tmp = add_te_fam(tmp, loop_te_data)
tmp = add_te_spec(tmp, enr_te_fams)
dat.K562 = tmp

# Combine the data for plotting
tmp.CH12 = table(c(dat.CH12$te_name_l, dat.CH12$te_name_r))
tmp.CH12 = tmp.CH12[order(tmp.CH12)]
tmp.CH12 = tmp.CH12[which(names(tmp.CH12) %in% enr_te_fams$name)]
tmp.GM12878 = table(c(dat.GM12878$te_name_l, dat.GM12878$te_name_r))
tmp.GM12878 = tmp.GM12878[order(tmp.GM12878)]
tmp.GM12878 = tmp.GM12878[which(names(tmp.GM12878) %in% enr_te_fams$name)]
tmp.K562 = table(c(dat.K562$te_name_l, dat.K562$te_name_r))
tmp.K562 = tmp.K562[order(tmp.K562)]
tmp.K562 = tmp.K562[which(names(tmp.K562) %in% enr_te_fams$name)]

te_names = union(union(names(tmp.K562), names(tmp.GM12878)), names(tmp.CH12))
dat = data.frame(name = te_names, species = rep("Non-Enriched", length(te_names)), CH12 = rep(NA, length(te_names)), GM12878 = rep(NA, length(te_names)), K562 = rep(NA, length(te_names)), stringsAsFactors=FALSE)
rownames(dat) = dat$name
dat[rownames(enr_te_fams),"species"] = enr_te_fams$species
dat[names(tmp.CH12),"CH12"] = tmp.CH12
dat[names(tmp.GM12878),"GM12878"] = tmp.GM12878
dat[names(tmp.K562),"K562"] = tmp.K562
library(reshape)
tmp = melt(dat)
tmp$name = factor(tmp$name, levels=dat$name)
tmp$species = factor(tmp$species, levels=c("Mouse", "Human", "Shared", "Non-Enriched"))




# Supplemental figure 3: Full TE overlap data for enriched families
pdf("loop_te-enriched-fam-intersections.pdf")
ggplot(tmp, aes(x=variable, y=value, fill=name)) +
geom_bar(stat="identity", position="fill")
dev.off()

# Figure 2D: TE Enrichment Types by Cell
pdf("loop_te-enriched-species-intersections.pdf")
ggplot(tmp, aes(x=variable, y=value, fill=species)) +
geom_bar(stat="identity", position="fill")
dev.off()
