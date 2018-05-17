all_ctcf_repeats$estAge = NA

# Substitution rates are from estimates presented in Kumar, S., & Subramanian, S. (2002). Mutation rates in mammalian genomes. Proc Natl Acad Sci U S A, 99(2), 803–808.
all_ctcf_repeats[which(all_ctcf_repeats$species == "hg19"),"estAge"] = (all_ctcf_repeats[which(all_ctcf_repeats$species == "hg19"),"pctDiv_rmsk"]/100) / 2.2e-9
all_ctcf_repeats[which(all_ctcf_repeats$species == "mm9"),"estAge"] = (all_ctcf_repeats[which(all_ctcf_repeats$species == "mm9"),"pctDiv_rmsk"]/100) / 2.4e-9

# Prepare the box plot of age distributions for binomial-enriched TEs.
dat = all_ctcf_repeats[which(all_ctcf_repeats$name_rmsk %in% rownames(dat)),c(1,5,15,20)]
row_order$V1 = rev(row_order$V1)
dat$order = factor(row_order[dat$name_rmsk,1])

pdf("Figure-1E_enriched-te-ages.boxplot.binomial.pdf")
ggplot(dat, aes(order, estAge/1000000)) +
geom_boxplot(aes(fill=species)) +
scale_x_discrete(labels=rownames(row_order), name="TE Type") +
scale_y_continuous(name="Millions of Years") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
dev.off()
