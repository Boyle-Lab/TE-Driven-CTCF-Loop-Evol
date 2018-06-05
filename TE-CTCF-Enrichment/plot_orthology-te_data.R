library(VennDiagram)

dat.href = read.table("CTCF.hg19.mm9.union.labelled.te.bed", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(dat.href) = c("chrom", "chromStart", "chromEnd", "ID", "peak", "Human", "Mouse", "TE")

pdf("Figure-1A_CTCF-occupancy_TE-overlap.venn.pdf")

draw.triple.venn(length(which(dat.href$Human == 1)), length(which(dat.href$Mouse == 1)), length(which(dat.href$TE == 1)), length(which(dat.href$Human == 1 & dat.href$Mouse == 1)), length(which(dat.href$Mouse == 1 & dat.href$TE == 1)), length(which(dat.href$Human == 1 & dat.href$TE == 1)), length(which(dat.href$Human == 1 & dat.href$Mouse == 1 & dat.href$TE == 1)), category=c("Human", "Mouse", "TE"), euler.d = TRUE, scaled = TRUE)

dev.off()
