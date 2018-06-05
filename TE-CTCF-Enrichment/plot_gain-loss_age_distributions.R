dat.te_gl = read.table("gain-loss_TE.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(dat.te_gl) = c("chrom", "start", "end", "id", "peak", "status", "TE", "name", "class", "family", "pctDiv", "species")
dat.te_gl[which(dat.te_gl$species == "hg19"),"estAge"] = (dat.te_gl[which(dat.te_gl$species == "hg19"),"pctDiv"]/100)/2.2e-9
dat.te_gl[which(dat.te_gl$species == "mm9"),"estAge"] = (dat.te_gl[which(dat.te_gl$species == "mm9"),"pctDiv"]/100)/2.4e-9
dat.te_gl[with(dat.te_gl, grepl("gain", status)),"gl"] = "gain"
dat.te_gl[with(dat.te_gl, grepl("loss", status)),"gl"] = "loss"
dat.te_gl[with(dat.te_gl, grepl("ortholog", status)),"gl"] = "ortholog"


library(ggplot2)

pdf("gain-loss_age-distr.pdf")
ggplot(dat.te_gl, aes(gl, estAge/1000000)) +
geom_boxplot(aes(fill=species), notch=TRUE) +
scale_y_continuous(name="Millions of Years") +
scale_x_discrete(name="Evolutinary History")
dev.off()