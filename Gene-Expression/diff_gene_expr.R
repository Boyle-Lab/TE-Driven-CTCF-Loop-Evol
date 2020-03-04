# R pipeline to produce histone mod and gene expression plots in Figure 7 and Sup. Fig. 13.

require("heatmap2")

# Read in data on loop annotations.
load("../TE-Loop-Intersection/loop_intersection.Rdata")

# Read in gene expression data.
gene_expr_data = read.table("../data/gene_exp/rna_seq/normCountsComplete.tab", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c(1,2,4:6)]
colnames(gene_expr_data) = c("gene_mm9", "gene_hg19", "CH12", "K562", "GM12878")

# Load up functions
source("add_bigwig_annotation.R")
source("annotate_enhancer_promoters.R")

# Annotate loop data with histone modifications at anchor loci.
tmp = add_all_histmods(annotated_loop_data, bigwig_path="../data/histmods", mods=c("H3K4me3", "H3K4me1", "H3K27me3", "H3K27ac"))
# Apply min-max normalization to histmod annotations
hist_annotated_loop_data = norm_all_histmods(tmp, mods=c("H3K4me3", "H3K4me1", "H3K27me3", "H3K27ac"))

# Produce histone mod plots for enhancer-enhancer, promoter-enhancer, and
# promoter-promoter interactions. (Sup. Fig. 13A)
dat = hist_annotated_loop_data[which(hist_annotated_loop_data$cell_q == "GM12878" & ((hist_annotated_loop_data$H3K4me3.max.l >= 5 & hist_annotated_loop_data$H3K4me1.max.r >= 7) | (hist_annotated_loop_data$H3K4me3.max.l >= 5 & hist_annotated_loop_data$H3K4me1.max.r >= 7) | (hist_annotated_loop_data$H3K4me1.max.l >= 7 & hist_annotated_loop_data$H3K4me1.max.r >= 7) | (hist_annotated_loop_data$H3K4me3.max.l >= 5 & hist_annotated_loop_data$H3K4me3.max.r >= 5)) ),]
for (i in 1:nrow(dat)) {  if (dat[i,47] > dat[i,46]) { dat[i,c(46,48,47,49)] = dat[i,c(47,49,46,48)] }}
pdf("eh-histmods.GM12878.pdf", height=10, width=5)
heatmap.2(as.matrix(dat[,c(46,48,47,49)]), trace="none", dendrogram="none", Colv=FALSE, col=brewer.pal(9,"Reds"), cexCol=1)
dev.off()

dat = hist_annotated_loop_data[which(hist_annotated_loop_data$cell_q == "K562" & ((hist_annotated_loop_data$H3K4me3.max.l >= 30 & hist_annotated_loop_data$H3K4me1.max.r >= 9) | (hist_annotated_loop_data$H3K4me3.max.l >= 30 & hist_annotated_loop_data$H3K4me1.max.r >= 9) | (hist_annotated_loop_data$H3K4me1.max.l >= 9 & hist_annotated_loop_data$H3K4me1.max.r >= 9) | (hist_annotated_loop_data$H3K4me3.max.l >= 30 & hist_annotated_loop_data$H3K4me3.max.r >= 30)) & hist_annotated_loop_data$H3K4me3.max.l <= 165 & hist_annotated_loop_data$H3K4me3.max.r <= 165),]
for (i in 1:nrow(dat)) {  if (dat[i,47] > dat[i,46]) { dat[i,c(46,48,47,49)] = dat[i,c(47,49,46,48)] }}
pdf("eh-histmods.K562.pdf", height=10, width=5)
heatmap.2(as.matrix(dat[,c(46,48,47,49)]), trace="none", dendrogram="none", Colv=FALSE, col=brewer.pal(9,"Reds"), cexCol=1)
dev.off()

dat = hist_annotated_loop_data[which(hist_annotated_loop_data$cell_q == "CH12" & ((hist_annotated_loop_data$H3K4me3.max.l >= 12 & hist_annotated_loop_data$H3K4me1.max.r >= 2) | (hist_annotated_loop_data$H3K4me3.max.l >= 12 & hist_annotated_loop_data$H3K4me1.max.r >= 2) | (hist_annotated_loop_data$H3K4me1.max.l >= 2 & hist_annotated_loop_data$H3K4me1.max.r >= 2) | (hist_annotated_loop_data$H3K4me3.max.l >= 12 & hist_annotated_loop_data$H3K4me3.max.r >= 12)) & hist_annotated_loop_data$H3K4me1.max.l <= 7 & hist_annotated_loop_data$H3K4me1.max.r <= 7  & hist_annotated_loop_data$H3K4me3.max.l <= 200 & hist_annotated_loop_data$H3K4me3.max.r <= 200),]
for (i in 1:nrow(dat)) {  if (dat[i,47] > dat[i,46]) { dat[i,c(46,48,47,49)] = dat[i,c(47,49,46,48)] }}
pdf("eh-histmods.CH12.pdf", height=10, width=5)
heatmap.2(as.matrix(dat[,c(46,48,47,49)]), trace="none", dendrogram="none", Colv=FALSE, col=brewer.pal(9,"Reds"), cexCol=1)
dev.off()

# Annotate loop anchors with their nearest TSS.
gene_hist_annotated_loop_data = annotate_anchor_tss(hist_annotated_loop_data)

# Test for differences in delta-TPM between conserved and nonconserved loops
# for each pair of cell types.

# GM12878 to K562
tmp.n = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "GM12878" & gene_hist_annotated_loop_data$class_t2 != "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
tmp.c = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "GM12878" & gene_hist_annotated_loop_data$class_t2 == "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
delta_exp.n = unlist(apply(tmp.n, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
delta_exp.c = unlist(apply(tmp.c, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
# Expected stats:
#nrow(tmp.n)
#[1] 2976
#nrow(tmp.c)
#[1] 1436
#
#summary(delta_exp.n)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.000    0.101    1.373   20.305   14.397 2055.215
#
#summary(delta_exp.c)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's
#   0.0000    0.0499    0.7218   14.5510    8.5328 2055.2150         1
#
#wilcox.test(delta_exp.n, delta_exp.c, alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  delta_exp.n and delta_exp.c
##W = 2325213, p-value = 8.055e-07
#alternative hypothesis: true location shift is greater than 0


# GM12878 to CH12
tmp.n = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "GM12878" & gene_hist_annotated_loop_data$class_t1 != "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
tmp.c = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "GM12878" & gene_hist_annotated_loop_data$class_t1 == "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
delta_exp.n = unlist(apply(tmp.n, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))
delta_exp.c = unlist(apply(tmp.c, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))
# Expected stats:
#length(delta_exp.c)
#[1] 63
#length(delta_exp.n)
#[1] 4349
#
#summary(delta_exp.n)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's
#   0.000    0.124    1.282   29.523   14.743 4084.146        1
#summary(delta_exp.c)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#  0.00000   0.04133   0.29874   8.26716   5.64290 134.48700
#
#wilcox.test(delta_exp.n, delta_exp.c, alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  delta_exp.n and delta_exp.c
#W = 160216, p-value = 0.01019
#alternative hypothesis: true location shift is greater than 0

# K562 to GM12878
tmp.n = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "K562" & gene_hist_annotated_loop_data$class_t2 != "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
tmp.c = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "K562" & gene_hist_annotated_loop_data$class_t2 == "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
delta_exp.n = unlist(apply(tmp.n, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
delta_exp.c = unlist(apply(tmp.c, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
# Expected stats:
#summary(delta_exp.n)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.0000    0.0699    1.0109   18.3965   10.8295 2055.2150
#summary(delta_exp.c)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's
#   0.0000    0.0510    0.7312   14.8593    8.6922 2055.2150         1
#length(delta_exp.c)
#[1] 1438
#length(delta_exp.n)
#[1] 1353
#wilcox.test(delta_exp.n, delta_exp.c, alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  delta_exp.n and delta_exp.c
#W = 1011256, p-value = 0.03279
#alternative hypothesis: true location shift is greater than 0

# K562 to CH12
tmp.n = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "K562" & gene_hist_annotated_loop_data$class_t1 != "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
tmp.c = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "K562" & gene_hist_annotated_loop_data$class_t1 == "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
delta_exp.n = unlist(apply(tmp.n, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))
delta_exp.c = unlist(apply(tmp.c, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))
# Expected stats:
#length(delta_exp.c)
#[1] 45
#length(delta_exp.n)
#[1] 2746
#summary(delta_exp.n)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's
#   0.000    0.082    1.205   25.942   15.457 3947.808        1
#summary(delta_exp.c)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#  0.00000   0.04127   0.32399  17.46553   7.09201 296.44700
#wilcox.test(delta_exp.n, delta_exp.c, alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  delta_exp.n and delta_exp.c
#W = 70671, p-value = 0.04806
#alternative hypothesis: true location shift is greater than 0

# CH12 to GM12878
tmp.n = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "CH12" & gene_hist_annotated_loop_data$class_t1 != "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
tmp.c = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "CH12" & gene_hist_annotated_loop_data$class_t1 == "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
delta_exp.n = unlist(apply(tmp.n, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))
delta_exp.c = unlist(apply(tmp.c, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))

# Expected stats:
#summary(delta_exp.n)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.0000   0.1969   1.3902  24.0994  13.6744 935.9437
#summary(delta_exp.c)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#  0.00000   0.05141   0.51529   8.36820   6.31739 134.48700
#wilcox.test(delta_exp.n, delta_exp.c, alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  delta_exp.n and delta_exp.c
#W = 7963, p-value = 0.01471
#alternative hypothesis: true location shift is greater than 0

# CH12 to K562
tmp.n = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "CH12" & gene_hist_annotated_loop_data$class_t2 != "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
tmp.c = gene_hist_annotated_loop_data[which(gene_hist_annotated_loop_data$cell_q == "CH12" & gene_hist_annotated_loop_data$class_t2 == "C" & ( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) )) ,]
delta_exp.n = unlist(apply(tmp.n, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
delta_exp.c = unlist(apply(tmp.c, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
# Expected stats:
#summary(delta_exp.n)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.0000   0.1587   2.1077  26.1118  22.0336 895.8701
#summary(delta_exp.c)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.0000   0.0445   0.3811  19.1285  11.1703 296.4470
#wilcox.test(delta_exp.n, delta_exp.c, alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  delta_exp.n and delta_exp.c
#W = 6743, p-value = 0.02512
#alternative hypothesis: true location shift is greater than 0

# Full delta-TPM comparison graphs for Sup. Figs 13B-C
dat = gene_hist_annotated_loop_data[which( (abs(gene_hist_annotated_loop_data$tss_dist.l) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.r) >= 5000) | (abs(gene_hist_annotated_loop_data$tss_dist.r) <= 1000 & abs(gene_hist_annotated_loop_data$tss_dist.l) >= 5000) ),]
dat$delta_exp.t1 = unlist(apply(dat, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t1.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t1.r"])) ) } } ))
dat$cons.t1 = unlist(apply(dat, 1, function(x) { if (x["class_t1"] == "C") { return("C") } else { return("N") } } ))
dat$delta_exp.t2 = unlist(apply(dat, 1, function(x) { if (x["tss_dist.l"] < x["tss_dist.r"]) { return( abs(as.numeric(x["expr_q.l"]) - as.numeric(x["expr_t2.l"])) ) } else { return( abs(as.numeric(x["expr_q.r"]) - as.numeric(x["expr_t2.r"])) ) } } ))
dat$cons.t2 = unlist(apply(dat, 1, function(x) { if (x["class_t2"] == "C") { return("C") } else { return("N") } } ))

# Compile into a long-format table
pdat = data.frame("cell_q"=dat$cell_q, "cell_t"=dat$cell_t1, "cons"=dat$cons.t1, "delta_exp"=dat$delta_exp.t1, "te_derived"=dat$te_derived, stringsAsFactors=FALSE)
pdat = rbind(pdat, data.frame("cell_q"=dat$cell_q, "cell_t"=dat$cell_t2, "cons"=dat$cons.t2, "delta_exp"=dat$delta_exp.t2, "te_derived"=dat$te_derived, stringsAsFactors=FALSE))
pdat$comp = paste(pdat$cell_q, pdat$cell_t, sep="-to-")

# Sup. Fig. 13B
pdf("delta_exp_cons-v-noncons.full.te.pdf")
x = ggplot(pdat[which(pdat$te_derived==TRUE),], aes(comp))
x = x + geom_bar(aes(y=delta_exp, fill=cons), position="dodge", stat="summary", fun.y="mean", na.rm=TRUE)
x
dev.off()

# Sup. Fig. 13C
pdf("delta_exp_cons-v-noncons.full.non-te.pdf")
x = ggplot(pdat[which(pdat$te_derived==FALSE),], aes(comp))
x = x + geom_bar(aes(y=delta_exp, fill=cons), position="dodge", stat="summary", fun.y="mean", na.rm=TRUE)
x
dev.off()

# Combined plot for Fig. 7A (Created in two parts)
pdat$sp_comp = apply(pdat, 1, function(x, species_idx=c("GM12878"="Human", "K562"="Human", "CH12"="Mouse")) { return(paste(species_idx[x["cell_q"]], species_idx[x["cell_t"]], sep="-")) } )

# Further simplify

pdat[which(pdat$sp_comp == "Human-Mouse"),"sp_comp"] = "Mouse-Human"

# TE-derived panel
pdf("delta_exp_cons-v-noncons.summary.te.pdf")
x = ggplot(pdat[which(pdat$te_derived==TRUE),], aes(sp_comp))
x = x + geom_bar(aes(y=delta_exp, fill=cons), position="dodge", stat="summary", fun.y="mean", na.rm=TRUE)
x + geom_point(data=pdat[which(pdat$te_derived==TRUE & pdat$cons=="C" & pdat$sp_comp=="Mouse-Human"),], aes(sp_comp, delta_exp, color=cons), position=position_jitter(width=0.02, height=0.02))
dev.off()

# Non-TE panel
pdf("delta_exp_cons-v-noncons.summary.non-te.pdf")
x = ggplot(pdat[which(pdat$te_derived==FALSE),], aes(sp_comp))
x = x + geom_bar(aes(y=delta_exp, fill=cons), position="dodge", stat="summary", fun.y="mean", na.rm=TRUE)
x
dev.off()

# Significance tests between sets with expected results and
# observation count.
# Human-Human Non-TE
#wilcox.test(pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Human-Human" & pdat$cons=="N"),"delta_exp"], pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Human-Human" & pdat$cons=="C"),"delta_exp"], alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  pdat[which(pdat$te_derived == FALSE & pdat$sp_comp == "Human-Human" &#  and pdat[which(pdat$te_derived == FALSE & pdat$sp_comp == "Human-Human" & #    pdat$cons == "N"), "delta_exp"] and     pdat$cons == "C"), "delta_exp"]
#W = 4351314, p-value = 2.849e-06
#alternative hypothesis: true location shift is greater than 0
#nrow(pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Human-Human" & pdat$cons=="N"),])
#[1] 3356
#nrow(pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Human-Human" & pdat$cons=="C"),])
#[1] 2424

# Human-Human TE-derived
#wilcox.test(pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Human-Human" & pdat$cons=="N"),"delta_exp"], pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Human-Human" & pdat$cons=="C"),"delta_exp"], alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  pdat[which(pdat$te_derived == TRUE & pdat$sp_comp == "Human-Human" & # and pdat[which(pdat$te_derived == TRUE & pdat$sp_comp == "Human-Human" &   #  pdat$cons == "N"), "delta_exp"] and     pdat$cons == "C"), "delta_exp"]
#W = 243500, p-value = 0.0001866
#alternative hypothesis: true location shift is greater than 0
#nrow(pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Human-Human" & pdat$cons=="C"),])
#[1] 450
#nrow(pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Human-Human" & pdat$cons=="N"),])
#[1] 973


# Mouse-Human Non-TE
#wilcox.test(pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="N"),"delta_exp"], pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="C"),"delta_exp"], alt="g")
#    Wilcoxon rank sum test with continuity correction
#data:  pdat[which(pdat$te_derived == FALSE & pdat$sp_comp == "Mouse-Human" &#  and pdat[which(pdat$te_derived == FALSE & pdat$sp_comp == "Mouse-Human" & #    pdat$cons == "N"), "delta_exp"] and     pdat$cons == "C"), "delta_exp"]
#W = 755044, p-value = 6.449e-05
#alternative hypothesis: true location shift is greater than 0
#nrow(pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="C"),])
#[1] 221
#nrow(pdat[which(pdat$te_derived==FALSE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="N"),])
#[1] 5935

# Mouse-Human TE-derived
#wilcox.test(pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="N"),"delta_exp"], pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="C"),"delta_exp"], alt="g")
#    Wilcoxon rank sum test with continuity correction]
#data:  pdat[which(pdat$te_derived == TRUE & pdat$sp_comp == "Mouse-Human" & # and pdat[which(pdat$te_derived == TRUE & pdat$sp_comp == "Mouse-Human" &   #  pdat$cons == "N"), "delta_exp"] and     pdat$cons == "C"), "delta_exp"]
#W = 8306.5, p-value = 0.1811
#alternative hypothesis: true location shift is greater than 0
nrow(pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="N"),])
#[1] 1572
nrow(pdat[which(pdat$te_derived==TRUE & pdat$sp_comp=="Mouse-Human" & pdat$cons=="C"),])
#[1] 9

# Heat map of CAMK2D expression values. (Fig. 7B)
pdf("gene_exp_vals.CAM2KD.pdf")
heatmap.2(as.matrix(cbind(rep(0,3),t(gene_expr_data[which(gene_expr_data$gene_hg19 == "CAMK2D"),c(3,5,4)]))),  trace="none", dendrogram="none", Colv=FALSE, col=colorRampPalette(c("white","red")), cexCol=1)
dev.off()
