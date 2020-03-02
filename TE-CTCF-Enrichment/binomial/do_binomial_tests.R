# Read in complete repeatMasker data
all_repeats = read.table("../data/hadoop/repeatmasker_with_tss_dists.dat", stringsAsFactors=FALSE, header=FALSE)
colnames(all_repeats) = c("id_rmsk", "chrom", "chromStart", "chromEnd", "name", "class", "family", "pctDiv", "species", "nearest_gene", "tss_dist")
all_repeats$tss_dist_bin = cut(all_repeats$tss_dist, breaks = c(-Inf,-10000,-2000,0,2000,10000,Inf), include.lowest = FALSE, labels=FALSE)
rownames(all_repeats) = all_repeats$id_rmsk

# Read in the CTCF-TE intersection data
all_ctcf_repeats = read.table("ctcf_te_intersection.all.txt", stringsAsFactors=FALSE, header=FALSE)
colnames(all_ctcf_repeats) = c("id_rmsk", "chrom", "chromStart_rmsk", "chromEnd_rmsk", "name_rmsk", "class_rmsk", "family_rmsk", "pctDiv_rmsk", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "factor", "species", "cell")

# Add TSS distance bins to CTCF-TE intersections
library(sqldf)
all_ctcf_repeats = sqldf("SELECT y.*, x.nearest_gene, x.tss_dist, x.tss_dist_bin FROM all_repeats x INNER JOIN all_ctcf_repeats y ON x.id_rmsk = y.id_rmsk")

# Calculate the expected genomic frequencies for each species/cell.
ctcf_exp_freqs = as.data.frame(matrix(nrow=1, ncol=4))
colnames(ctcf_exp_freqs) = c("GM12878", "K562", "CH12", "MEL")
for (cell in c("GM12878", "K562", "CH12", "MEL")) {
    species = "hg19"
    if (cell == "CH12" || cell == "MEL") {
        species = "mm9"
    }
    ctcf_exp_freqs[,cell] = nrow(all_ctcf_repeats[which(all_ctcf_repeats$cell == cell),]) / nrow(all_repeats[which(all_repeats$species == species),])
}

# Perform binomial enrichment tests for each TE family with CTCF overlaps.
ctcf_te_enr = as.data.frame(matrix(nrow=length(unique(all_ctcf_repeats$name_rmsk)), ncol=24))
colnames(ctcf_te_enr) = c("te_count.GM12878", "te_ctcf_count.GM12878", "obs_freq.GM12878", "exp_freq.GM12878", "pbinom.GM12878", "pbinom_cor.GM12878", "te_count.K562", "te_ctcf_count.K562", "obs_freq.K562", "exp_freq.K562", "pbinom.K562", "pbinom_cor.K562", "te_count.CH12", "te_ctcf_count.CH12", "obs_freq.CH12", "exp_freq.CH12", "pbinom.CH12", "pbinom_cor.CH12", "te_count.MEL", "te_ctcf_count.MEL", "obs_freq.MEL", "exp_freq.MEL", "pbinom.MEL", "pbinom_cor.MEL")
rownames(ctcf_te_enr) = unique(all_ctcf_repeats$name_rmsk)

for (cell in c("GM12878", "K562", "CH12", "MEL")) {

    species = "hg19"
    if (cell == "CH12" || cell == "MEL") {
        species = "mm9"
    }

    te_c_col = paste("te_count", cell, sep=".")
    te_ctcf_col = paste("te_ctcf_count", cell, sep=".")
    of_col = paste("obs_freq", cell, sep=".")
    ef_col = paste("exp_freq", cell, sep=".")
    p_col = paste("pbinom", cell, sep=".")
    pc_col = paste("pbinom_cor", cell, sep=".")
    ctcf_te_enr[, ef_col] = ctcf_exp_freqs[,cell]

    for (rp in unique(all_ctcf_repeats$name_rmsk)) {
        ctcf_te_enr[rp, te_c_col] = nrow(all_repeats[which(all_repeats$name == rp & all_repeats$species == species),])
        ctcf_te_enr[rp, te_ctcf_col] = nrow(all_ctcf_repeats[which(all_ctcf_repeats$name == rp & all_ctcf_repeats$cell == cell),])
        ctcf_te_enr[rp, of_col] = ctcf_te_enr[rp, te_ctcf_col] / ctcf_te_enr[rp, te_c_col]
        if (is.na(ctcf_te_enr[rp, te_c_col]) || ctcf_te_enr[rp, te_c_col] == 0) {
            #print(rp)
            next
        }

        f = binom.test(ctcf_te_enr[rp, te_ctcf_col], ctcf_te_enr[rp, te_c_col], alternative = "g", p = ctcf_exp_freqs[,cell])
        ctcf_te_enr[rp, p_col] = f$p.value
    }
    ctcf_te_enr[which(!is.na(ctcf_te_enr[,te_c_col])), pc_col] = p.adjust(ctcf_te_enr[which(!is.na(ctcf_te_enr[,te_c_col])), p_col], method="holm")
}

####
# Supplemental table 2: complete binomial enrichments
# Write all rows where at least one cell meets the following criteria:
# 1) p-value <= 1e-4
# 2) >= 25 CTCF-bound TE insertions
write.table(ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),], quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", file="Supplemental-Table-2_ctcf_te_enrichments.binomial.txt")


####
# Figure 2b: Heat map of binomial enrichments.
library(ggplot2)
library(reshape)
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
dat = 1-dat
dat$family = rownames(dat)
dat $order = factor(dat$family, levels=(dat$family)[order(dat$pbinom_cor.GM12878, dat$pbinom_cor.K562, dat$pbinom_cor.CH12, dat$pbinom_cor.MEL)])
dat = melt(dat)
dat$value[which(is.na(dat$value))] = 0
pdf("Fig-2b_te-enrichment-heatmap.binomial.pdf")
ggplot(dat, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()
# Note that some reordering of rows was performed by hand in the presentation figure.
# Final row ordering vector is stored in row-order.txt.

####
# Figure 2c: Heat map of CTCF-bound insertion count for enriched families.
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(2,8,14,20)]
dat = log2(dat)
library(ggplot2)
library(reshape)
dat$family = rownames(dat)
row_order = read.table("row-order.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
rownames(row_order) = row_order[,2]
dat$order = factor(row_order[rownames(dat),1])
dat.melt = melt(dat)
dat.melt$value[which(is.na(dat.melt$value))] = 0
pdf("Fig-2c_te-count-heatmap.binomial.pdf")
ggplot(dat.melt, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()


####
# Supplemental figure 1c: Euler diagram of enriched families by cell.
library(venneuler)
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24)]
dat[is.na(dat)] = 1
dat[dat > 0.05] = 1
dat[dat<1]=0
dat = abs(dat-1)
vd = venneuler(dat)
pdf("Sup-fig-1c_te-enriched-families_euler.binomial.pdf")
plot(vd)
dev.off()


####
# Figure 2d: binding frequencies for all enriched TE families.
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(3,9,15,21)]
dat = -log2(dat)
dat$family = rownames(dat)
dat$order = factor(row_order[rownames(dat),1])
dat.melt = melt(dat)
dat.melt$value[which(is.na(dat.melt$value))] = 0
pdf("Fig-2d_te-frequencies-heatmap.binomial.pdf")
ggplot(dat.melt, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()

####
# Figure 2e: Species-specific enrichment pie charts.

# enriched_te_fams.txt was hand-curated to annotate species-specificity
# of TE type enrichment.
enr_te_fams = read.table("enriched_te_fams.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, col.names=c("species", "name"), row.names="name")

all_ctcf_repeats$te_spec = enr_te_fams[all_ctcf_repeats$name_rmsk,"species"]
all_ctcf_repeats[which(is.na(all_ctcf_repeats$te_spec)),"te_spec"] = "Non-Enriched"
dat = melt(table(all_ctcf_repeats[,c("species","te_spec")]))
dat$te_spec = factor(dat$te_spec, levels=c("Mouse", "Shared", "Human", "Non-Enriched"))

pdf("frac_ctcf_enr.pie.hg19.pdf")
ggplot(dat[which(dat$species=="hg19"),], aes(x=species, y=value, fill=te_spec)) +
geom_bar(stat="identity", position="fill") +
coord_polar("y", start=0)
dev.off()

pdf("frac_ctcf_enr.pie.mm9.pdf")
ggplot(dat[which(dat$species=="mm9"),], aes(x=species, y=value, fill=te_spec)) +
geom_bar(stat="identity", position="fill") +
coord_polar("y", start=0)
dev.off()

####
# Figure 2f: Motif scores in enriched, non-enriched, random repeats, and background seqs
all_motif_scores = read.table("all-repeats_CTCF-ren.motif_scores.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)

# motif scores for enriched and non-enriched repeats
all_ancestral_motif_scores = all_motif_scores[which(all_motif_scores$V1 %in% unique(all_ctcf_repeats$name_rmsk)),]
all_ancestral_motif_scores$max = unlist(lapply(all_ancestral_motif_scores[,6], function(str) { x = as.numeric(unlist(strsplit(str, ","))); return(max(x, na.rm=TRUE))}))
all_ancestral_motif_scores$source = "consensus"
all_ancestral_motif_scores[which(!(all_ancestral_motif_scores$V1 %in% enr_te_fams$name)),"source"] = "Non-Enriched"
all_ancestral_motif_scores[which(all_ancestral_motif_scores$V1 %in% enr_te_fams$name),"source"] = "Enriched"

# Background sequences
library("rtracklayer")
library("regioneR")
all_ancestral_motif_scores$length = unlist(lapply(all_ancestral_motif_scores[,6], function(str) { return(length(as.numeric(unlist(strsplit(str, ",")))))}))
bg_regions = createRandomRegions(nregions=nrow(all_ancestral_motif_scores), length.mean = mean(all_ancestral_motif_scores$length), length.sd = sd(all_ancestral_motif_scores$length), genome = "hg19", non.overlapping = TRUE)
tmp = import("../..//data/motifs/hg19.CTCF.ren.bigWig", selection=bg_regions, as="NumericList")
tmp = unlist(lapply(tmp, max))
dat = all_ancestral_motif_scores[,c("source", "max")]
dat = rbind(dat, data.frame("source" = rep("Background", length(tmp)), "max" = tmp))

# Random repeats
rand_motif_scores = all_repeat_scores[sample(all_repeat_scores$V1, nrow(all_ancestral_motif_scores)),]
rand_motif_scores$max = unlist(lapply(rand_motif_scores[,6], function(str) { x = as.numeric(unlist(strsplit(str, ","))); return(max(x, na.rm=TRUE))}))
rand_motif_scores$source = "random_repeats"
dat = rbind(dat, rand_motif_scores[,c("source", "max")])
dat$source = factor(dat$source, levels=c("Enriched", "Non-Enriched", "random_repeats", "Background"))

# Plot the data
pdf("CTCF-bound-repeat-consensus_motif-scores.pdf")
ggplot(dat, aes(x=source, y=max, fill=source)) +
geom_boxplot(notch=TRUE)
dev.off()

# Save data structures used in other analyses.
save(all_repeats, all_ctcf_repeats, ctcf_te_enr, row_order, file="binomial_enrichments.Rdata")
