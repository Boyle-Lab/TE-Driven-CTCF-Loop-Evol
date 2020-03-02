# Plot CTCF motif score distributions for HOA repeats (Fig. 3C)

require("rtracklayer")
require("reshape")

# Load enrichment test data.
load("../binomial/binomial_enrichments.Rdata")

# Load the gain-loss data.
dat.te_gl = read.table("gain-loss_TE.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(dat.te_gl) = c("chrom", "start", "end", "id", "peak", "status", "TE", "name", "class", "family", "pctDiv", "species")
dat.te_gl[with(dat.te_gl, grepl("gain", status)),"gl"] = "gain"
dat.te_gl[with(dat.te_gl, grepl("loss", status)),"gl"] = "loss"
dat.te_gl[with(dat.te_gl, grepl("ortholog", status)),"gl"] = "ortholog"

# Build a Genomic Ranges object of 50bp windows surrounding CTCF peaks embedded in mouse repeats.
ints = data.frame(chrom=all_ctcf_repeats$chrom, start=all_ctcf_repeats$peak_ctcf-25, end=all_ctcf_repeats$peak_ctcf+24, id=all_ctcf_repeats$id_ctcf, cell=all_ctcf_repeats$cell)
ints = ints[which(ints$cell %in% c("CH12", "MEL")),]
ints$chrom = as.character(ints$chrom)
ints = GRanges(ints, mcols=data.frame(id=ints[,"id"]))

# Select a matching number of 50bp background regions, excluding the bound regions...
bg_regions = createRandomRegions(nregions=length(ints), length.mean=50, length.sd=0, genome = "mm9", non.overlapping = TRUE, mask = ints)

# Get the maximum motif scores within these regions.
tmp = import("../../data/motifs/mm9.CTCF.ren.bigWig", selection=bg_regions, as="NumericList")
tmp = unlist(lapply(tmp, max))

# Get the human regions as well.
ints = data.frame(chrom=all_ctcf_repeats$chrom, start=all_ctcf_repeats$peak_ctcf-25, end=all_ctcf_repeats$peak_ctcf+24, id=all_ctcf_repeats$id_ctcf, cell=all_ctcf_repeats$cell)
ints = ints[which(ints$cell %in% c("K562", "GM12878")),]
ints$chrom = as.character(ints$chrom)
ints = GRanges(ints, mcols=data.frame(id=ints[,"id"]))

# Select a matching number of 50bp background regions, excluding the bound regions.
bg_regions = createRandomRegions(nregions=length(ints), length.mean=50, length.sd=0, genome = "hg19", non.overlapping = TRUE, mask = ints)

# Get the maximum CTCF motif scores.
tmp1 = import("../../data/motifs/hg19.CTCF.ren.bigWig", selection=bg_regions, as="NumericList")
tmp1 = unlist(lapply(tmp, max))

# Assemble plot data.
dat = melt(dat.te_gl[which(dat.te_gl$status == "ortholog" & dat.te_gl$name %in% human_only),c("name", "mscore_hg19", "mscore_mm9")])

dat$name = "Human-Only"
dat = rbind(dat, data.frame(name = rep("Human-Only",nrow(human_only_motif_scores)), variable = rep("consensus",nrow(human_only_motif_scores)), value = human_only_motif_scores$max))

#  We need 852 total random regions to match the observed data. Select an equal number of human and mouse regions.
tmp2 = sample(tmp, 852/2)
tmp2 = append(tmp2, sample(tmp1, 852/2))

# Add these to the plot data.
dat = rbind(dat, data.frame(name = rep("Human-Only",length(tmp2)), variable = rep("Background",length(tmp2)), value = tmp2))
dat$variable = as.character(dat$variable)
dat[which(dat$variable == "mscore_hg19"),"variable"] = "Human"
dat[which(dat$variable == "mscore_mm9"),"variable"] = "Mouse"
dat[which(dat$variable == "consensus"),"variable"] = "Consensus"
dat$variable = factor(dat$variable, levels=c("Human", "Mouse", "Consensus", "Background"))

# Draw the box plots.
pdf("enr_te_motif_scores.human-only_with-cons-and-rand.pdf")
ggplot(dat, aes(x=name, y=value)) +
geom_boxplot(aes(fill=variable),notch=TRUE) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()