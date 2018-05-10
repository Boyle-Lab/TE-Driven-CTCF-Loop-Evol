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
# Figure 1b: Heat map of binomial enrichments.
library(ggplot2)
library(reshape)
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
dat = 1-dat
dat$family = rownames(dat)
dat $order = factor(dat$family, levels=(dat$family)[order(dat$pbinom_cor.GM12878, dat$pbinom_cor.K562, dat$pbinom_cor.CH12, dat$pbinom_cor.MEL)])
dat = melt(dat)
dat$value[which(is.na(dat$value))] = 0
pdf("Fig-1a_te-enrichment-heatmap.binomial.pdf")
ggplot(dat, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()
# Note that some reordering of rows was performed by hand in the presentation figure.
# Final row ordering vector is stored in row-order.txt.

####
# Figure 1c: Heat map of CTCF-bound insertion count for enriched families.
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
pdf("Fig-1b_te-count-heatmap.binomial.pdf")
ggplot(dat.melt, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()


####
# Figure 1d: Euler diagram of enriched families by cell.
library(venneuler)
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24)]
dat[is.na(dat)] = 1
dat[dat > 0.05] = 1
dat[dat<1]=0
dat = abs(dat-1)
vd = venneuler(dat)
pdf("Fig-1c_te-enriched-families_euler.binomial.pdf")
plot(vd)
dev.off()


####
# Supplemental figure 1: binding frequencies for all enriched TE families.
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(3,9,15,21)]
dat = -log2(dat)
dat$family = rownames(dat)
dat$order = factor(row_order[rownames(dat),1])
dat.melt = melt(dat)
dat.melt$value[which(is.na(dat.melt$value))] = 0
pdf("Supplemental-fig-1_te-frequencies-heatmap.binomial.pdf")
ggplot(dat.melt, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()

####
# Table 1: Top-5 enriched TE families for shared, human, and mouse

# Shared
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) & (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) & (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 10 & ctcf_te_enr$obs_freq.CH12 >= 0.01) & (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 10 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24,2,8,14,20,3,9,15,21)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
binomial_top5 = dat[1:5,]

# Human-specific
dat = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) & (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01)  & (is.na(ctcf_te_enr$pbinom_cor.CH12)  | ctcf_te_enr$pbinom_cor.CH12 > 0.05) & (is.na(ctcf_te_enr$pbinom_cor.MEL)  | ctcf_te_enr$pbinom_cor.MEL > 0.05)),c(6,12,18,24,2,8,14,20,3,9,15,21)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
binomial_top5 = rbind(binomial_top5, dat[1:5,])

# Mouse-specific
dat = ctcf_te_enr[which((is.na(ctcf_te_enr$pbinom_cor.GM12878) | ctcf_te_enr$pbinom_cor.GM12878 > 0.05) & (is.na(ctcf_te_enr$pbinom_cor.K562) | ctcf_te_enr$pbinom_cor.K562 > 0.05) & (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) &(ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24,2,8,14,20,3,9,15,21)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
binomial_top5 = rbind(binomial_top5, dat[1:5,])
write.table(binomial_top5, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", file="Table1_top-5-ctcf-te-enrichments-per-species.txt")


