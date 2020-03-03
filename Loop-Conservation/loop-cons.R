# Assembles data from loop cross-mapping into data tables.

require("ggplot2")
require("reshape")
require("VennDiagram")

## Load up te-loop intersection data
load("../TE-Loop-Intersection/loop_intersection.Rdata")

## Load in commands
source("intersect_cons_te.R")

## Read in data

# Loop conservation classification results
loop_conservation_data = read.table("RAD21_loops.all.dat", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names(loop_conservation_data) = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "IAB", "FDR", "mapped_chr_l", "mapped_start_l", "mapped_end_l", "qstrand_l", "mapped_chr_r", "mapped_start_r", "mapped_end_r", "qstrand_r", "id_l", "anchor_l", "chrom_l", "start_l", "end_l",  "id_r", "anchor_r", "chrom_r", "start_r", "end_r", "class", "cell_q", "cell_t"))

## Intersect the conservation data with the TE-loop intersections
dat.GM12878 = intersect_cons_te(loop_conservation_data, dat.GM12878, "CH12", "K562")
dat.K562 = intersect_cons_te(loop_conservation_data, dat.K562, "CH12", "GM12878")
dat.CH12 = intersect_cons_te(loop_conservation_data, dat.CH12, "GM12878", "K562")

## Compile tabular data to incorporate into supplementary table 6
# First put all cell-wise data into a single list...  
dat = mget(ls(pattern = "dat.[A-Z]"))
names(dat) = lapply(names(dat), function (e) {unlist(strsplit(e, '.', fixed=TRUE))[2]}) 
write.table(compile_summary_data(dat), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, file="TE-derived-loop-cons.txt")


## Generate stacked bar charts for supplementary figure 6
dat.classes = read.table("RAD21_loops.stats.dat", header=TRUE, stringsAsFactors=FALSE, sep = "\t")
dat.classes$B2 = dat.classes$B2S+dat.classes$B2D
dat.classes = dat.classes[,c(1,5,10,7,6,3,4,2)]

pdf("class-breakdown-by-cell-comparison.scaled.pdf")
ggplot(melt(dat.classes), aes(x=File, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()


## Generate the species-comparison-wise stacked bar charts for Figure 5B
tmp = data.frame(matrix(unlist(c((dat.classes[1,2:8] + dat.classes[2,2:8]), (dat.classes[3,2:8] + dat.classes[5,2:8]), (dat.classes[4,2:8] + dat.classes[6,2:8])), nrow=3, byrow=T))
colnames(tmp) = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")
tmp$Comparison =  c("mouse-human", "human-mouse", "human-human")

pdf("class-breakdown-by-species-comparison.scaled.pdf")
ggplot(melt(tmp), aes(x=Comparison, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()


## Generate the species-wise bar charts for TAD overlaps (Fig. 5C).
TAD_data = read.table("../data/TAD/rao-and-huntley-2014_arrowhead-domains.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE, col.names=c("chrom", "x", "y", "score", "species", "cell"))
# Use 10kb resolution to call TAD overlaps and leave out human-mouse comparison since
# the mouse Hi-C data are much less complete than human datasets.
hits = get_tad_overlaps(loop_conservation_data[which(loop_conservation_data$cell_t != "CH12"),], TAD_data, max_dist = 10000)
# Rejigger the data for plotting.
x = table(hits[which(hits$is_tad == TRUE),c("cell_q", "class")])
x = melt(x)
x$class = factor(x$class, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
# Draw the plot
dat = rbind("Mouse" = x["CH12",], "Human" = x["GM12878",] + x["K562",])
dat = melt(dat)
colnames(dat) = c("species", "class", "value")
dat$class = factor(dat$class, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("TAD-conservation_by-class.pdf")
ggplot(dat, aes(x=species, y=value, fill=class)) +
geom_bar(stat="identity", position="fill")
dev.off()


## Supplementary Fig. 7 show conservation class membership with
## varying loop conservation resolution.
# 10KB
dat.classes.10k = read.table("../loop_conservation/RAD21_loops.stats.10k.dat", header=TRUE, stringsAsFactors=FALSE, sep = "\t")
dat.classes.10k$B2 = dat.classes.10k$B2S+dat.classes.10k$B2D
dat.classes.10k = dat.classes.10k[,c(1,5,10,7,6,3,4,2)]
tmp = data.frame(matrix( unlist( c( (dat.classes.10k[1,2:8] + dat.classes.10k[2,2:8]), (dat.classes.10k[3,2:8] + dat.classes.10k[5,2:8]), (dat.classes.10k[4,2:8] + dat.classes.10k[6,2:8]))), nrow=3, byrow=T))
colnames(tmp) = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")
tmp$Comparison =  factor(c("mouse-human", "human-mouse", "human-human"), levels=c("mouse-human", "human-mouse", "human-human"))
pdf("class-breakdown-by-species-comparison.10k.scaled.pdf")
ggplot(melt(tmp), aes(x=Comparison, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()

# 20KB
dat.classes.20k = read.table("../loop_conservation/RAD21_loops.stats.20k.dat", header=TRUE, stringsAsFactors=FALSE, sep = "\t")
dat.classes.20k$B2 = dat.classes.20k$B2S+dat.classes.20k$B2D
dat.classes.20k = dat.classes.20k[,c(1,5,10,7,6,3,4,2)]
tmp = data.frame(matrix( unlist( c( (dat.classes.20k[1,2:8] + dat.classes.20k[2,2:8]), (dat.classes.20k[3,2:8] + dat.classes.20k[5,2:8]), (dat.classes.20k[4,2:8] + dat.classes.20k[6,2:8]))), nrow=3, byrow=T))
colnames(tmp) = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")
tmp$Comparison =  factor(c("mouse-human", "human-mouse", "human-human"), levels=c("mouse-human", "human-mouse", "human-human"))
pdf("class-breakdown-by-species-comparison.20k.scaled.pdf")
ggplot(melt(tmp), aes(x=Comparison, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()

# 50KB
dat.classes.50k = read.table("../loop_conservation/RAD21_loops.stats.50k.dat", header=TRUE, stringsAsFactors=FALSE, sep = "\t")
dat.classes.50k = read.table("RAD21_loops.stats.50k.dat", header=TRUE, stringsAsFactors=FALSE, sep = "\t")
dat.classes.50k$B2 = dat.classes.50k$B2S+dat.classes.50k$B2D
dat.classes.50k = dat.classes.50k[,c(1,5,10,7,6,3,4,2)]
tmp = data.frame(matrix( unlist( c( (dat.classes.50k[1,2:8] + dat.classes.50k[2,2:8]), (dat.classes.50k[3,2:8] + dat.classes.50k[5,2:8]), (dat.classes.50k[4,2:8] + dat.classes.50k[6,2:8]))), nrow=3, byrow=T))
colnames(tmp) = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")
tmp$Comparison =  factor(c("mouse-human", "human-mouse", "human-human"), levels=c("mouse-human", "human-mouse", "human-human"))
pdf("class-breakdown-by-species-comparison.50k.scaled.pdf")
ggplot(melt(tmp), aes(x=Comparison, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()


# Produce the Venn diagrams in Supplementary Fig. 8
pdf("Loop-TAD-overlap.Human.pdf")
draw.pairwise.venn(length(which(hits$cell_q != "CH12")), length(which(TAD_data$cell != "CH12")), length(which(hits$is_tad == TRUE & hits$cell_q != "CH12")), category=c("Loops", "TADs"), scaled=TRUE)
dev.off()

pdf("Loop-TAD-overlap.Mouse.pdf")
draw.pairwise.venn(length(which(hits$cell_q == "CH12")), length(which(TAD_data$cell == "CH12")), length(which(hits$is_tad == TRUE & hits$cell_q == "CH12")), category=c("Loops", "TADs"), scaled=TRUE)
dev.off()

pdf("TAD-Conservation_Dixon-metric.Human-Human.pdf")
draw.pairwise.venn(length(which(hits$is_tad == TRUE & hits$cell_q == "GM12878" & hits$cell_t == "K562")), length(which(hits$is_tad == TRUE & hits$cell_q == "K562" & hits$cell_t == "GM12878")), length(which(hits$is_tad == TRUE & hits$cell_q == "K562" & hits$cell_t == "GM12878" & hits$class %in% c("C", "B2", "B1", "N1A"))), category=c("GM12878 TADs", "K562 TADs"), scaled=TRUE)
dev.off()

pdf("TAD-Conservation_Dixon-metric.Mouse-Human.pdf")
draw.pairwise.venn(length(which(hits$is_tad == TRUE & hits$cell_q == "CH12")), length(which(hits$is_tad == TRUE)), length(which(hits$is_tad == TRUE & hits$cell_q == "CH12" & hits$class %in% c("C", "B2", "B1", "N1A"))), category=c("Mouse TADs", "Human TADs"), scaled=TRUE)
dev.off()


## Figure 5D shows the spatial distribution of phastCons conservation scores across 500 bp
## windows surrounding the annotated CTCF ChIP-seq peak in TE-derived loop anchors, broken
## down by conservation class.
loop_ctcf_data = read.table("../data/hadoop/results/ctcf_loop_intersection.hg19.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("id_cp", "chrom", "start1", "end1", "start2", "end2", "iab", "fdr", "cell", "species", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak")
test2 = compile_class_bigwig_scores(dat.GM12878, loop_ctcf_data, "../data/phastCons/hg19/phastCons46way.placental.bigWig", "CH12", size=500)
# Convert data to long format matrix. For some reason, the "melt" function doesn't play well here...
dat2 = convert_to_long(test2)
# Plot the data as lines...
pdf("loop-conservation_vs_phastcons.GM12878-to-CH12.pdf")
ggplot(dat2, aes(x=position, y=value, colour=cat))+
geom_line()
dev.off()

# Mouse-human shows the same trends, but noisier due to the small set sizes.
test = compile_class_bigwig_scores(dat.CH12, hic_ctcf_data, "../data/phastCons/mm9/phastCons30way.placental.bigWig", "GM12878", "K562", size=1000)
dat = convert_to_long(test)
ggplot(dat, aes(x=position, y=value, colour=cat))+
geom_line()

# Human-human data have similar distributions for all classes, with only slight conservation drop at the peaks.
test3 = compile_class_bigwig_scores(dat.GM12878, loop_ctcf_data, "../data/phastCons/hg19/phastCons46way.placental.bigWig", "K562", size=1000, cats=c("C","B2","B1","B0"))
dat3 = convert_to_long(test3, cats=c("C","B2","B1","B0"))
ggplot(dat3, aes(x=position, y=value, colour=cat))+
geom_line()


## Supplemental Figure 9 contains the PhastCons score plots for TE-Derived loops only,
## showing the same score trends are evident.
test2 = compile_class_bigwig_scores(dat.GM12878, loop_te_data, "../data/phastCons/hg19/phastCons46way.placental.bigWig", "CH12", size=500)
dat2 = convert_to_long(test2)
# Plot the data as lines...
pdf("Supplemental_Figure-8_loop-conservation-wrt-phastCons_TE-derived-loops-only.pdf")
ggplot(dat2, aes(x=position, y=value, colour=cat))+
geom_line()
dev.off()

