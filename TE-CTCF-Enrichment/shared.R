# Read in complete repeatMasker data
all_repeats = read.table("../data/hadoop/repeatmasker_with_tss_dists.dat", stringsAsFactors=FALSE, header=FALSE)
colnames(all_repeats) = c("id_rmsk", "chrom", "chromStart", "chromEnd", "name", "class", "family", "pctDiv", "species", "nearest_gene", "tss_dist")
all_repeats$tss_dist_bin = cut(all_repeats$tss_dist, breaks = c(-Inf,-10000,-2000,0,2000,10000,Inf), include.lowest = FALSE, labels=FALSE)
rownames(all_repeats) = all_repeats$id_rmsk

# Read in the CTCF-TE intersection data
all_ctcf_repeats = read.table("ctcf_te_intersection.all.txt", stringsAsFactors=FALSE, header=FALSE)
colnames(all_ctcf_repeats) = c("id_rmsk", "chrom", "chromStart_rmsk", "chromEnd_rmsk", "name_rmsk", "class_rmsk", "family_rmsk", "pctDiv_rmsk", "id_ctcf", "chromStart_ctcf", "chromEnd_ctcf", "signalValue_ctcf", "peak_ctcf", "factor", "species", "cell")

# Add TSS distance bins to CTCF-TE intersections
library(sqldf)
all_ctcf_repeats = sqldf("SELECT y.*, x.nearest_gene, x.tss_dist, x.tss_dist_bin FROM all_repeats x INNER JOIN all_ctcf_repeats y ON x.id_rmsk = y.id_rmsk")
