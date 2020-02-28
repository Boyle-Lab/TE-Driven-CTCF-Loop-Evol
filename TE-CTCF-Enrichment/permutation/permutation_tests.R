# Do the permutation-based TE enrichment tests.
source("do_shuffle_peaks.R")
ctcf_te_enr_perm = do_permutations(all_ctcf_repeats, all_repeats)

# Assemble the data into a wide-format table.
ctcf_te_enr_perm_wide = as.data.frame(matrix(nrow=length(unique(ctcf_te_enr_perm$name)), ncol=32))
rownames(ctcf_te_enr_perm_wide) = unique(ctcf_te_enr_perm$name)
colnames(ctcf_te_enr_perm_wide)[1:8] = paste(colnames(ctcf_te_enr_perm), "GM12878", sep=".")[1:8]
colnames(ctcf_te_enr_perm_wide)[9:16] = paste(colnames(ctcf_te_enr_perm), "K562", sep=".")[1:8]
colnames(ctcf_te_enr_perm_wide)[17:24] = paste(colnames(ctcf_te_enr_perm), "CH12", sep=".")[1:8]
colnames(ctcf_te_enr_perm_wide)[25:32] = paste(colnames(ctcf_te_enr_perm), "MEL", sep=".")[1:8]
ctcf_te_enr_perm_wide[unique(ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "GM12878"), "name"]),1:8] = ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "GM12878"),c(1:8)]
ctcf_te_enr_perm_wide[unique(ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "K562"), "name"]),9:16] = ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "K562"),c(1:8)]
ctcf_te_enr_perm_wide[unique(ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "CH12"), "name"]),17:24] = ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "CH12"),c(1:8)]
ctcf_te_enr_perm_wide[unique(ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "MEL"), "name"]),25:32] = ctcf_te_enr_perm[which(ctcf_te_enr_perm$cell == "MEL"),c(1:8)]

# Add observed frequencies
ctcf_te_enr_perm_wide$obs_freq.GM12878 = ctcf_te_enr[rownames(ctcf_te_enr_perm_wide),"obs_freq.GM12878"]
ctcf_te_enr_perm_wide$obs_freq.K562 = ctcf_te_enr[rownames(ctcf_te_enr_perm_wide),"obs_freq.K562"]
ctcf_te_enr_perm_wide$obs_freq.CH12 = ctcf_te_enr[rownames(ctcf_te_enr_perm_wide),"obs_freq.CH12"]
ctcf_te_enr_perm_wide$obs_freq.MEL = ctcf_te_enr[rownames(ctcf_te_enr_perm_wide),"obs_freq.MEL"]


# Supplemental table 3: Full permuation enrichment data.
# Criteria for enrichments:
# 1) p-value <= 1e-4
# 2) >= 25 CTCF-bound insertions
# 3) CTCF binding frequency >= 0.01
dat = ctcf_te_enr_perm_wide[which((ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.GM12878 >= 25 & ctcf_te_enr_perm_wide$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.K562 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.K562 >= 25 & ctcf_te_enr_perm_wide$obs_freq.K562 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.CH12 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.CH12 >= 10 & ctcf_te_enr_perm_wide$obs_freq.CH12 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.MEL <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.MEL >= 10 & ctcf_te_enr_perm_wide$obs_freq.MEL >= 0.01)),]
dat = cbind(dat, ctcf_te_enr[rownames(dat),c(1,7,13,19)])
write.table(dat, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", file="ctcf_te_enrichments.permutation.strict.txt")


# Supplemental figure 2a: Enrichment p-value heat map
dat = ctcf_te_enr_perm_wide[which((ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.GM12878 >= 25 & ctcf_te_enr_perm_wide$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.K562 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.K562 >= 25 & ctcf_te_enr_perm_wide$obs_freq.K562 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.CH12 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.CH12 >= 10 & ctcf_te_enr_perm_wide$obs_freq.CH12 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.MEL <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.MEL >= 10 & ctcf_te_enr_perm_wide$obs_freq.MEL >= 0.01)),c(8,16,24,32)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
dat = 1-dat
library(ggplot2)
library(reshape)
dat$family = rownames(dat)
dat $order = factor(dat$family, levels=(dat$family)[order(dat$pval_bn_cor.GM12878, dat$pval_bn_cor.K562, dat$pval_bn_cor.CH12, dat$pval_bn_cor.MEL)])
dat.melt = melt(dat)
dat.melt$value[which(is.na(dat.melt$value))] = 0
pdf("Supplemental-Figure-2a_TE-repeat-enrichments_permuation-pvals.pdf")
ggplot(dat.melt, aes(x=variable, y=order, fill=value)) + geom_tile() + scale_fill_gradient(high=rgb(0.44,0.1,0.62), low="white") + theme(axis.text.y = element_text(size=4))
dev.off()


# Supplemental figure 2b: Euler diagram of enriched TE families by cell.
library(venneuler)
dat = ctcf_te_enr_perm_wide[which((ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.GM12878 >= 25 & ctcf_te_enr_perm_wide$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.K562 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.K562 >= 25 & ctcf_te_enr_perm_wide$obs_freq.K562 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.CH12 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.CH12 >= 10 & ctcf_te_enr_perm_wide$obs_freq.CH12 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.MEL <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.MEL >= 10 & ctcf_te_enr_perm_wide$obs_freq.MEL >= 0.01)),c(8,16,24,32)]
dat[is.na(dat)] = 1
dat[dat > 0.05] = 1
dat[dat<1]=0
dat = abs(dat-1)
vd = venneuler(dat)
pdf("Supplemental-Figure-2b_enriched-te-families-by-cell.permuation.pdf")
plot(vd)
dev.off()


# Supplemental figure 2c: Euler diagram showing overlap between binomial and permutation-based enrichment tests.
dat1 = ctcf_te_enr_perm_wide[which((ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.GM12878 >= 25 & ctcf_te_enr_perm_wide$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.K562 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.K562 >= 25 & ctcf_te_enr_perm_wide$obs_freq.K562 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.CH12 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.CH12 >= 10 & ctcf_te_enr_perm_wide$obs_freq.CH12 >= 0.01) | (ctcf_te_enr_perm_wide$pval_bn_cor.MEL <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.MEL >= 10 & ctcf_te_enr_perm_wide$obs_freq.MEL >= 0.01)),c(8,16,24,32)]
dat2 = ctcf_te_enr[which((ctcf_te_enr$pbinom_cor.GM12878 <= 0.0001 & ctcf_te_enr$te_ctcf_count.GM12878 >= 25 & ctcf_te_enr$obs_freq.GM12878 >= 0.01) | (ctcf_te_enr$pbinom_cor.K562 <= 0.0001 & ctcf_te_enr$te_ctcf_count.K562 >= 25 & ctcf_te_enr$obs_freq.K562 >= 0.01) | (ctcf_te_enr$pbinom_cor.CH12 <= 0.0001 & ctcf_te_enr$te_ctcf_count.CH12 >= 25 & ctcf_te_enr$obs_freq.CH12 >= 0.01) | (ctcf_te_enr$pbinom_cor.MEL <= 0.0001 & ctcf_te_enr$te_ctcf_count.MEL >= 25 & ctcf_te_enr$obs_freq.MEL >= 0.01)),c(6,12,18,24)]
tmp = union(rownames(dat1), rownames(dat2))
dat = as.data.frame(matrix(nrow=length(tmp), ncol=2))
rownames(dat) = tmp
colnames(dat) = c("Binomial", "Permutation")
at[] = 0
dat[rownames(dat2),1] = 1
dat[rownames(dat1),2] = 1
vd = venneuler(dat)
pdf("Supplemental-figure-2c_enriched-te-families-permutation-vs-binomial.pdf")
plot(vd)
dev.off()


# Supplemental figure 2d: Top-5 enriched TE families for permuation-based tests

# Shared
dat = ctcf_te_enr_perm_wide[which((ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.GM12878 >= 25 & ctcf_te_enr_perm_wide$obs_freq.GM12878 >= 0.01) & (ctcf_te_enr_perm_wide$pval_bn_cor.K562 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.K562 >= 25 & ctcf_te_enr_perm_wide$obs_freq.K562 >= 0.01) & (ctcf_te_enr_perm_wide$pval_bn_cor.CH12 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.CH12 >= 10 & ctcf_te_enr_perm_wide$obs_freq.CH12 >= 0.01) & (ctcf_te_enr_perm_wide$pval_bn_cor.MEL <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.MEL >= 10 & ctcf_te_enr_perm_wide$obs_freq.MEL >= 0.01)),c(8,16,24,32,1,9,17,25,33,34,35,36)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
perm_top5 = dat[1:5,]

# Human-specific
dat = ctcf_te_enr_perm_wide[which((ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.GM12878 >= 25 & ctcf_te_enr_perm_wide$obs_freq.GM12878 >= 0.01) & (ctcf_te_enr_perm_wide$pval_bn_cor.K562 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.K562 >= 25 & ctcf_te_enr_perm_wide$obs_freq.K562 >= 0.01) & (is.na(ctcf_te_enr_perm_wide$pval_bn_cor.CH12) | ctcf_te_enr_perm_wide$pval_bn_cor.CH12 > 0.05) & (is.na(ctcf_te_enr_perm_wide$pval_bn_cor.MEL) | ctcf_te_enr_perm_wide$pval_bn_cor.CH12 > 0.05)),c(8,16,24,32,1,9,17,25,33,34,35,36)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
perm_top5 = rbind(perm_top5, dat[1:5,])

# Mouse-specific
dat = ctcf_te_enr_perm_wide[which((is.na(ctcf_te_enr_perm_wide$pval_bn_cor.GM12878) | ctcf_te_enr_perm_wide$pval_bn_cor.GM12878 > 0.05)  & (is.na(ctcf_te_enr_perm_wide$pval_bn_cor.K562) | ctcf_te_enr_perm_wide$pval_bn_cor.K562 > 0.05) & (ctcf_te_enr_perm_wide$pval_bn_cor.CH12 <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.CH12 >= 25 & ctcf_te_enr_perm_wide$obs_freq.CH12 >= 0.01) & (ctcf_te_enr_perm_wide$pval_bn_cor.MEL <= 0.0001 & ctcf_te_enr_perm_wide$observed_count.MEL >= 25 & ctcf_te_enr_perm_wide$obs_freq.MEL >= 0.01)),c(8,16,24,32,1,9,17,25,33,34,35,36)]
dat = dat[order(dat[,1], dat[,2], dat[,3], dat[,4]),]
perm_top5 = rbind(perm_top5, dat[1:5,])
write.table(perm_top5, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", file="Supplemental-figure-2d_ctcf-te-enrichments.permutation.top5.txt")



