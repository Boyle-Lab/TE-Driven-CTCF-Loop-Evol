# Calculate odds ratios for species-specificity of CTCF binding words.

# Load the whole-genome CTCF motif predictions
CTCF_genomic.hg19 = read.table("../../data/motifs/FIMO/CTCF.hg19.bed", header=FALSE, stringsAsFactors=FALSE, sep="\t")
CTCF_genomic.mm9 = read.table("../../data/motifs/FIMO/CTCF.mm9.bed", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(CTCF_genomic.hg19) = c("chrom", "start", "end", "name", "score", "strand", "word")
colnames(CTCF_genomic.mm9) = c("chrom", "start", "end", "name", "score", "strand", "word")

# Load the CTCF predictions for bound repeats
bound_rmsk_motifs.hg19 = read.table("CTCF.bound_rmsk.hg19.bed", header=FALSE, stringsAsFactors=FALSE, sep="\t")
bound_rmsk_motifs.mm9 = read.table("CTCF.bound_rmsk.mm9.bed", header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(bound_rmsk_motifs.hg19) = c("rpt", "start", "stop", "name", "score", "strand", "word")
colnames(bound_rmsk_motifs.mm9) = c("rpt", "start", "stop", "name", "score", "strand", "word")
bound_rmsk_motifs.hg19$word = factor(bound_rmsk_motifs.hg19$word)
bound_rmsk_motifs.mm9$word = factor(bound_rmsk_motifs.mm9$word)

# Summarize the motif-word content for human and mouse individually, and for the dataset as a whole
word_counts_rmsk.hg19 = summary(bound_rmsk_motifs.hg19$word, maxsum = nrow(bound_rmsk_motifs.hg19))
word_counts_rmsk.mm9 = summary(bound_rmsk_motifs.mm9$word, maxsum = nrow(bound_rmsk_motifs.mm9))
bound_rmsk_motifs.all = rbind(bound_rmsk_motifs.hg19, bound_rmsk_motifs.mm9)
word_counts_rmsk.all.df = as.data.frame(summary(bound_rmsk_motifs.all$word, maxsum = nrow(bound_rmsk_motifs.all)))
colnames(word_counts_rmsk.all.df)[1] = "occ.all"
word_counts_rmsk.all.df[,"occ.hg19"] = word_counts_rmsk.hg19[rownames(word_counts_rmsk.all.df)]
word_counts_rmsk.all.df[,"occ.mm9"] = word_counts_rmsk.mm9[rownames(word_counts_rmsk.all.df)]

# Calculate the normalization factors for each species, according to Schmidt et al.
# 20 is the length of the CTCF motif used in our analysis.
norm_factor.mm9 = (sum(word_counts_rmsk.mm9)*20)/1000000
norm_factor.hg19 = (sum(word_counts_rmsk.hg19)*20)/1000000

# Normalize counts
word_counts_rmsk.all.df$nocc.hg19 = word_counts_rmsk.all.df$occ.hg19 / norm_factor.hg19
word_counts_rmsk.all.df$nocc.mm9 = word_counts_rmsk.all.df$occ.mm9 / norm_factor.mm9
word_counts_rmsk.all.df[is.na(word_counts_rmsk.all.df)] = 0  # Zero out all missing values

# Calculate odds ratios for species-specificity. Positive values are human-specific and negative values are mouse-specific.
word_counts_rmsk.all.df$orr = log(word_counts_rmsk.all.df$nocc.hg19 / word_counts_rmsk.all.df$nocc.mm9)

# Calculate enrichment p-values for each family with individual fisher's exact tests, with bonferroni p-value adjustment at the end.
source("rpt_word_enrichments.R")
bound_rpt_word_enr.mm9 = calc_rpt_enrichments(bound_rmsk_motifs.mm9, word_counts_rmsk.all.df, CTCF_genomic.mm9, "mm9")
bound_rpt_word_enr.hg19 = calc_rpt_enrichments(bound_rmsk_motifs.hg19, word_counts_rmsk.all.df, CTCF_genomic.hg19, "hg19")

# Write the enriched TEs to Supplementary Table 4: species-specific-word-enriched-TEs.
dat = bound_rpt_word_enr.hg19[which(bound_rpt_word_enr.hg19$specific_fg > 0 & bound_rpt_word_enr.hg19$pval_cor <= 10e-40),]
dat = rbind(dat, bound_rpt_word_enr.mm9[which(bound_rpt_word_enr.mm9$specific_fg > 0 & bound_rpt_word_enr.mm9$pval_cor <= 10e-40),])
dat$species = c(rep("human", 5), rep("mouse", 5))
write.table(dat, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE, file="Supplemental-Table-4_species-specific-word-enriched-TEs.txt")

# For each of the TE types observed in our binomial tests, get the word-enrichment stats for the most frequent CTCF motif-word. This is presented in Supplementary Table 5.
source("get_top_enr.R")
binomial_enriched_te = read.table("row-order.txt", stringsAsFactors=FALSE, header=FALSE, sep="\t")
binomial_enriched_te = binomial_enriched_te$V2
tmp.hg = get_te_word_data(bound_rmsk_motifs.hg19, word_counts_rmsk.all.df, binomial_enriched_te, "hg19")
tmp.mm = get_te_word_data(bound_rmsk_motifs.mm9, word_counts_rmsk.all.df, binomial_enriched_te, "mm9")
write.table(cbind(tmp.hg, tmp.mm), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE, file="Supplemental-Table-5_top-enriched-words.binomial.txt")
