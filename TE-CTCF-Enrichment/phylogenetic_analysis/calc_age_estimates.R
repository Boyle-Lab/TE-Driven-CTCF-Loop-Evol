# Calculate age estimates for binomial TE enrichment dataset
# and produce plot for Figure 3A.

# Load up binomial enrichment dataset. This contains
# all_repeats, all_ctcf_repeats, ctcf_te_enr, and row_order objects
load("../binomial/binomial_enrichments.Rdata")

# Add age estimates based on pctDiv estimates
all_ctcf_repeats$estAge = NA
all_ctcf_repeats[which(all_ctcf_repeats$species == "hg19"),"estAge"] = (all_ctcf_repeats[which(all_ctcf_repeats$species == "hg19"),"pctDiv_rmsk"]/100) / 2.2e-9
all_ctcf_repeats[which(all_ctcf_repeats$species == "mm9"),"estAge"] = (all_ctcf_repeats[which(all_ctcf_repeats$species == "mm9"),"pctDiv_rmsk"]/100) / 2.4e-9

# Draw the age distribution box-plots used in Fig. 3A
dat = all_ctcf_repeats[which(all_ctcf_repeats$name_rmsk %in% rownames(dat)),c(1,5,15,20)]
row_order$V1 = rev(row_order$V1)
dat$order = factor(row_order[dat$name_rmsk,1])

pdf("enriched-te-ages.boxplot.binomial.pdf")
ggplot(dat, aes(order, estAge/1000000)) +
geom_boxplot(aes(fill=species)) +
scale_x_discrete(labels=rownames(row_order), name="TE Type") +
scale_y_continuous(name="Millions of Years") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
dev.off()

