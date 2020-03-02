# Produce plots for phylogenetic gain and loss or TE-associated CTCF binding sites.
require(ggplot2)
require(reshape)

# Load enrichment test data.
load("../binomial/binomial_enrichments.Rdata")

# Read in data from the gain-loss prediction pipeline.
dat.te_gl = read.table("gain-loss_TE.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(dat.te_gl) = c("chrom", "start", "end", "id", "peak", "status", "TE", "name", "class", "family", "pctDiv", "species")
dat.te_gl[with(dat.te_gl, grepl("gain", status)),"gl"] = "gain"
dat.te_gl[with(dat.te_gl, grepl("loss", status)),"gl"] = "loss"
dat.te_gl[with(dat.te_gl, grepl("ortholog", status)),"gl"] = "ortholog"

# Draw the bar charts in Figure 3B.
dat = table(dat.te_gl[which(dat.te_gl$name %in% enr_te_fams$name),c("name", "status")])
tmp = melt(dat)
tmp$status = factor(tmp$status, levels=c("ortholog", "loss_mm9", "gain_hg19", "loss_hg19", "gain_mm9"))
tmp$name = factor(tmp$name, levels=c("MamGypLTR3a", "LTR16", "X12_DNA", "X1_DNA", "LTR41C", "L1PA17", "L1PA15-16", "L1P4a", "HERVIP10F-int", "MER52A", "HUERS-P3-int", "HUERS-P3b-int", "HUERS-P2-int", "MER66C", "LTR8B", "HERVK-int", "MER1B", "LTR1A1", "L1PBa", "HERVIP10FH-int", "MER1A", "LTR1A2", "LTR1", "LTR13_", "MER61-int", "THE1B-int", "HERVE-int", "MER41A", "HERVL-int", "LTR13", "MER91A", "LTR16E1", "LTR16C", "LTR33A", "LTR33", "MER102c", "MER102b", "LTR33A_", "ERV3-16A3_I-int", "LTR16A1", "Tigger16a", "LTR16A2", "AmnSINE1", "MER20B", "MER102a", "MER74A", "MER119", "MER135", "L1M1", "L1M3f", "ERVL-B4-int", "MER52-int", "MamGypLTR3", "MER91B", "LTR41B", "LTR55", "LTR41", "LTR50", "MER91C", "MER20", "ORR1E", "B3A", "ORR1D2", "B3", "RMER1A", "RMER1B", "B2_Mm2","B2_Mm1t", "B2_Mm1a", "IAPEY4_LTR"))

pdf("enr-te-fams_gain-loss-orth.pdf")
ggplot(tmp, aes(x=name, y=value, fill=status)) +
geom_bar(stat="identity", position="fill") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()
