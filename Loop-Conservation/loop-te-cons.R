# R analysis pipeline for data presented in Fig. 6, Sup. Figs. 11-12, and Sup. Table 6.

require("ggplot2")
require("reshape")
require("UpSetR")

## Load up annotations from elsewhere.
load("../TE-Loop-Intersection/loop_intersection.Rdata")
load("loop-cons.Rdata")

## Load in commands
source("intersect_cons_te.R")
source("do_upset.R")

## Read in data

# Loop conservation classification results
loop_conservation_data = read.table("RAD21_loops.all.dat", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names(loop_conservation_data) = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "IAB", "FDR", "mapped_chr_l", "mapped_start_l", "mapped_end_l", "qstrand_l", "mapped_chr_r", "mapped_start_r", "mapped_end_r", "qstrand_r", "id_l", "anchor_l", "chrom_l", "start_l", "end_l",  "id_r", "anchor_r", "chrom_r", "start_r", "end_r", "class", "cell_q", "cell_t"))

## Compile tabular data to incorporate into supplementary table 6
# First put all cell-wise data into a single list...  
dat = mget(ls(pattern = "dat.[A-Z]"))
names(dat) = lapply(names(dat), function (e) {unlist(strsplit(e, '.', fixed=TRUE))[2]}) 
write.table(compile_summary_data(dat), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, file="TE-derived-loop-cons.txt")

## Figures 6A-C, and sup. fig 11 are hybrid Upset/Bar plots.
# This requires two parts to be drawn separately and combined:

# First prepare the set data.
dat.sets.kg = list("Conserved" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "C"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "C" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B2" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "B2"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "B2" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B1" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "B1"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "B1" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B0" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "B0"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "B0" & (dat.K562$te_left | dat.K562$te_right)),])),
                "N1A" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "N1A"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "N1A" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "N1B"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "N1B" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "N0"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "N0" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),]))
               )

dat.sets.kc = list("Conserved" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "C"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "C" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B2" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "B2"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "B2" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B1" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "B1"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "B1" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B0" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "B0"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "B0" & (dat.K562$te_left | dat.K562$te_right)),])),
                "N1A" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "N1A"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "N1A" & ( (dat.K562$te_left & dat.K562$CH12_l_maps == FALSE) | (dat.K562$te_right & dat.K562$CH12_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "N1B"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "N1B" & ( (dat.K562$te_left & dat.K562$CH12_l_maps == FALSE) | (dat.K562$te_right & dat.K562$CH12_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.K562[which(dat.K562$class_CH12 == "N0"),]), nrow(dat.K562[which(dat.K562$class_CH12 == "N0" & ( (dat.K562$te_left & dat.K562$CH12_l_maps == FALSE) | (dat.K562$te_right & dat.K562$CH12_r_maps == FALSE) )),]))
               )

dat.sets.gk = list("Conserved" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "C"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "C" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B2" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B2"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B2" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B1" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B1"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B1" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B0" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B0"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B0" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "N1A" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1A"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1A" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1B"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1B" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N0"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N0" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),]))
               )

dat.sets.gc = list("Conserved" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "C"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "C" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B2" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B2"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B2" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B1" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B1"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B1" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B0" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B0"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B0" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "N1A" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1A"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1A" & ( (dat.GM12878$te_left & dat.GM12878$CH12_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$CH12_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1B"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1B" & ( (dat.GM12878$te_left & dat.GM12878$CH12_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$CH12_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N0"),]), nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N0" & ( (dat.GM12878$te_left & dat.GM12878$CH12_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$CH12_r_maps == FALSE) )),]))
               )

dat.sets.ck = list("Conserved" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "C"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "C" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B2" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "B2"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "B2" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B1" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "B1"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "B1" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B0" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "B0"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "B0" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "N1A" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "N1A"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "N1A" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "N1B"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "N1B" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "N0"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "N0" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),]))
               )

dat.sets.cg = list("Conserved" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B2" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B2"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B2" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B1" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B1"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B1" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B0" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B0"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B0" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "N1A" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1A"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1A" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1B"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1B" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),]))
               )

# Generate the bar plots for Figs 6A-C
dat.sets.mm_hg = as.data.frame(t(as.data.frame(dat.sets.cg))) + as.data.frame(t(as.data.frame(dat.sets.ck)))
dat.sets.mm_hg[,1] = dat.sets.mm_hg[,1] - dat.sets.mm_hg[,2]
colnames(dat.sets.mm_hg) = c("No TE", "TE")
dat.sets.mm_hg$cat = factor(rownames(dat.sets.mm_hg), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_mouse-to-human.pdf")
ggplot(melt(dat.sets.mm_hg), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.hg_mm = as.data.frame(t(as.data.frame(dat.sets.gc))) + as.data.frame(t(as.data.frame(dat.sets.kc)))
dat.sets.hg_mm[,1] = dat.sets.hg_mm[,1] - dat.sets.hg_mm[,2]
colnames(dat.sets.hg_mm) = c("No TE", "TE")
dat.sets.hg_mm$cat = factor(rownames(dat.sets.hg_mm), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_human-to-mouse.pdf")
ggplot(melt(dat.sets.hg_mm), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.hg_hg = as.data.frame(t(as.data.frame(dat.sets.gk))) + as.data.frame(t(as.data.frame(dat.sets.kg)))
dat.sets.hg_hg[,1] = dat.sets.hg_hg[,1] - dat.sets.hg_hg[,2]
colnames(dat.sets.hg_hg) = c("No TE", "TE")
dat.sets.hg_hg$cat = factor(rownames(dat.sets.hg_hg), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_human-to-human.pdf")
ggplot(melt(dat.sets.hg_hg), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

# UpSet plots for Figs 6A-C
dat_upset.ck = c("Conserved" = nrow(dat.CH12[which(dat.CH12$class_K562 == "C"),]),
              "Conserved&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "C" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "B2" = nrow(dat.CH12[which(dat.CH12$class_K562 == "B2"),]),
              "B2&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "B2" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "B1" = nrow(dat.CH12[which(dat.CH12$class_K562 == "B1"),]),
              "B1&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "B1" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "B0" = nrow(dat.CH12[which(dat.CH12$class_K562 == "B0"),]),
              "B0&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "B0" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "N1A" = nrow(dat.CH12[which(dat.CH12$class_K562 == "N1A"),]),
              "N1A&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "N1A" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),]),
              "N1B" = nrow(dat.CH12[which(dat.CH12$class_K562 == "N1B"),]),
              "N1B&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "N1B" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),]),
              "N0" = nrow(dat.CH12[which(dat.CH12$class_K562 == "N0"),]),
              "N0&TE" = nrow(dat.CH12[which(dat.CH12$class_K562 == "N0" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),]))

dat_upset.cg = c("Conserved" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C"),]),
              "Conserved&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "B2" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B2"),]),
              "B2&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B2" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "B1" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B1"),]),
              "B1&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B1" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "B0" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B0"),]),
              "B0&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B0" & (dat.CH12$te_left | dat.CH12$te_right)),]),
              "N1A" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1A"),]),
              "N1A&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1A" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),]),
              "N1B" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1B"),]),
              "N1B&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1B" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),]),
              "N0" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0"),]),
              "N0&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),]))

dat_upset.kg = c("Conserved" = nrow(dat.K562[which(dat.K562$class_GM12878 == "C"),]),
              "Conserved&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "C" & (dat.K562$te_left | dat.K562$te_right)),]),
              "B2" = nrow(dat.K562[which(dat.K562$class_GM12878 == "B2"),]),
              "B2&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "B2" & (dat.K562$te_left | dat.K562$te_right)),]),
              "B1" = nrow(dat.K562[which(dat.K562$class_GM12878 == "B1"),]),
              "B1&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "B1" & (dat.K562$te_left | dat.K562$te_right)),]),
              "B0" = nrow(dat.K562[which(dat.K562$class_GM12878 == "B0"),]),
              "B0&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "B0" & (dat.K562$te_left | dat.K562$te_right)),]),
              "N1A" = nrow(dat.K562[which(dat.K562$class_GM12878 == "N1A"),]),
              "N1A&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "N1A" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),]),
              "N1B" = nrow(dat.K562[which(dat.K562$class_GM12878 == "N1B"),]),
              "N1B&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "N1B" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),]),
              "N0" = nrow(dat.K562[which(dat.K562$class_GM12878 == "N0"),]),
              "N0&TE" = nrow(dat.K562[which(dat.K562$class_GM12878 == "N0" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),]))

dat_upset.kc = c("Conserved" = nrow(dat.K562[which(dat.K562$class_CH12 == "C"),]),
              "Conserved&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "C" & (dat.K562$te_left | dat.K562$te_right)),]),
              "B2" = nrow(dat.K562[which(dat.K562$class_CH12 == "B2"),]),
              "B2&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "B2" & (dat.K562$te_left | dat.K562$te_right)),]),
              "B1" = nrow(dat.K562[which(dat.K562$class_CH12 == "B1"),]),
              "B1&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "B1" & (dat.K562$te_left | dat.K562$te_right)),]),
              "B0" = nrow(dat.K562[which(dat.K562$class_CH12 == "B0"),]),
              "B0&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "B0" & (dat.K562$te_left | dat.K562$te_right)),]), "N1A" = nrow(dat.K562[which(dat.K562$class_CH12 == "N1A"),]),
              "N1A&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "N1A" & ( (dat.K562$te_left & dat.K562$CH12_l_maps == FALSE) | (dat.K562$te_right & dat.K562$CH12_r_maps == FALSE) )),]),
              "N1B" = nrow(dat.K562[which(dat.K562$class_CH12 == "N1B"),]),
              "N1B&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "N1B" & ( (dat.K562$te_left & dat.K562$CH12_l_maps == FALSE) | (dat.K562$te_right & dat.K562$CH12_r_maps == FALSE) )),]),
              "N0" = nrow(dat.K562[which(dat.K562$class_CH12 == "N0"),]),
              "N0&TE" = nrow(dat.K562[which(dat.K562$class_CH12 == "N0" & ( (dat.K562$te_left & dat.K562$CH12_l_maps == FALSE) | (dat.K562$te_right & dat.K562$CH12_r_maps == FALSE) )),]))

dat_upset.gc = c("Conserved" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "C"),]),
              "Conserved&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "C" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "B2" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B2"),]),
              "B2&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B2" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "B1" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B1"),]),
              "B1&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B1" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "B0" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B0"),]),
              "B0&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "B0" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "N1A" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1A"),]),
              "N1A&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1A" & ( (dat.GM12878$te_left & dat.GM12878$CH12_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$CH12_r_maps == FALSE) )),]),
              "N1B" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1B"),]),
              "N1B&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N1B" & ( (dat.GM12878$te_left & dat.GM12878$CH12_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$CH12_r_maps == FALSE) )),]),
              "N0" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N0"),]),
              "N0&TE" = nrow(dat.GM12878[which(dat.GM12878$class_CH12 == "N0" & ( (dat.GM12878$te_left & dat.GM12878$CH12_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$CH12_r_maps == FALSE) )),]))

dat_upset.gk = c("Conserved" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "C"),]),
              "Conserved&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "C" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "B2" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B2"),]),
              "B2&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B2" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "B1" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B1"),]),
              "B1&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B1" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "B0" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B0"),]),
              "B0&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B0" & (dat.GM12878$te_left | dat.GM12878$te_right)),]),
              "N1A" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1A"),]),
              "N1A&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1A" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),]),
              "N1B" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1B"),]),
              "N1B&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1B" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),]),
              "N0" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N0"),]),
              "N0&TE" = nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N0" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),]))

dat_upset.hg_mm = dat_upset.gc + dat_upset.kc
dat_upset.hg_hg = dat_upset.gk + dat_upset.kg
dat_upset.mm_hg = dat_upset.cg + dat_upset.ck

pdf("upset_mouse-to-human.pdf")
upset(fromExpression(dat_upset.mm_hg), nsets=length(dat_upset.mm_hg), order.by="degree", group.by="sets")
dev.off()

pdf("upset_human-to-mouse.pdf")
upset(fromExpression(dat_upset.hg_mm), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()

pdf("upset_human-to-human.pdf")
upset(fromExpression(dat_upset.hg_hg), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()



# Bar plots for Sup. Fig. 11.
dat.sets.cg = as.data.frame(t(as.data.frame(dat.sets.cg)))
dat.sets.cg[,1] = dat.sets.cg[,1] - dat.sets.cg[,2]
colnames(dat.sets.cg) = c("No TE", "TE")
dat.sets.cg$cat = factor(rownames(dat.sets.cg), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_CH12-to-GM12878.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.ck = as.data.frame(t(as.data.frame(dat.sets.ck)))
dat.sets.ck[,1] = dat.sets.ck[,1] - dat.sets.ck[,2]
colnames(dat.sets.ck) = c("No TE", "TE")
dat.sets.ck$cat = factor(rownames(dat.sets.ck), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_CH12-to-K562.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.gk = as.data.frame(t(as.data.frame(dat.sets.gk)))
dat.sets.gk[,1] = dat.sets.gk[,1] - dat.sets.gk[,2]
colnames(dat.sets.gk) = c("No TE", "TE")
dat.sets.gk$cat = factor(rownames(dat.sets.gk), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_GM12878-to-K562.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.kg = as.data.frame(t(as.data.frame(dat.sets.kg)))
dat.sets.kg[,1] = dat.sets.kg[,1] - dat.sets.kg[,2]
colnames(dat.sets.kg) = c("No TE", "TE")
dat.sets.kg$cat = factor(rownames(dat.sets.kg), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_K562-to-GM12878.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.kc = as.data.frame(t(as.data.frame(dat.sets.kc)))
dat.sets.kc[,1] = dat.sets.kc[,1] - dat.sets.kc[,2]
colnames(dat.sets.kc) = c("No TE", "TE")
dat.sets.kc$cat = factor(rownames(dat.sets.kc), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_k562-to-ch12.pdf")
ggplot(melt(dat.sets.kc), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

dat.sets.gc = as.data.frame(t(as.data.frame(dat.sets.gc)))
dat.sets.gc[,1] = dat.sets.gc[,1] - dat.sets.gc[,2]
colnames(dat.sets.gc) = c("No TE", "TE")
dat.sets.gc$cat = factor(rownames(dat.sets.gc), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
pdf("stacked-bar_gm12878-to-ch12.pdf")
ggplot(melt(dat.sets.gc), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

# Generate the upset plot components...
pdf("upset_CH12-to-GM12878.pdf")
upset(fromExpression(dat_upset.cg), nsets=length(dat_upset.cg), order.by="degree", group.by="sets")
dev.off()

pdf("upset_CH12-to-K562.pdf")
upset(fromExpression(dat_upset.ck), nsets=length(dat_upset.ck), order.by="degree", group.by="sets")
dev.off()

pdf("upset_K562-to-GM12878.pdf")
upset(fromExpression(dat_upset.kg), nsets=length(dat_upset.kg), order.by="degree", group.by="sets")
dev.off()

pdf("upset_GM12878-to-K562.pdf")
upset(fromExpression(dat_upset.gk), nsets=length(dat_upset.gk), order.by="degree", group.by="sets")
dev.off()

pdf("upset_k562-to-ch12.pdf")
upset(fromExpression(dat_upset.kc), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()

pdf("upset_gm12878-to-ch12.pdf")
upset(fromExpression(dat_upset.gc), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()


## Sup. Fig. 12 examines resolution of conservation calls on TE contributions.

# 10k
loop_conservation_data.10k = read.table("../loop_conservation/RAD21_loops.all.10k.dat", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names=c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "IAB", "FDR", "mapped_chr_l", "mapped_start_l", "mapped_end_l", "qstrand_l", "mapped_chr_r", "mapped_start_r", "mapped_end_r", "qstrand_r", "id_l", "anchor_l", "chrom_l", "start_l", "end_l",  "id_r", "anchor_r", "chrom_r", "start_r", "end_r", "class", "cell_q", "cell_t"))
annotated_loop_data.10k = cbind(loop_data[which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"),c(1:4,6:7)], cell=rep("GM12878", length(which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"))))
annotated_loop_data.10k = rbind(annotated_loop_data.10k, cbind(loop_data[which(loop_data$cell == "K562" & loop_data$factor == "RAD21"),c(1:4,6:7)], cell=rep("K562", length(which(loop_data$cell == "K562" & loop_data$factor == "RAD21")))))
annotated_loop_data.10k = rbind(annotated_loop_data.10k, cbind(loop_data[which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C"),c(1:4,6:7)], cell=rep("CH12", length(which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C")))))
annotated_loop_data.10k = add_te_intersections(annotated_loop_data.10k, loop_te_data)
annotated_loop_data.10k = add_motif_orientations(annotated_loop_data.10k, loop_motif_data)
annotated_loop_data.10k = intersect_cons_te_2(loop_conservation_data.10k, annotated_loop_data.10k)
dat_upset =  do_upset_plot(annotated_loop_data.s25k, "CH12", "GM12878", "10k")
pdf("upset_CH12-to-GM12878.10k.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()
dat.sets =  do_bar_plot(annotated_loop_data.s25k, "CH12", "GM12878", "10k")
pdf("stacked-bar_CH12-to-GM12878.10k.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

# 20k
loop_conservation_data.20k = read.table("../loop_conservation/RAD21_loops.all.20k.dat", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names=c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "IAB", "FDR", "mapped_chr_l", "mapped_start_l", "mapped_end_l", "qstrand_l", "mapped_chr_r", "mapped_start_r", "mapped_end_r", "qstrand_r", "id_l", "anchor_l", "chrom_l", "start_l", "end_l",  "id_r", "anchor_r", "chrom_r", "start_r", "end_r", "class", "cell_q", "cell_t"))
annotated_loop_data.20k = cbind(loop_data[which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"),c(1:4,6:7)], cell=rep("GM12878", length(which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"))))
annotated_loop_data.20k = rbind(annotated_loop_data.20k, cbind(loop_data[which(loop_data$cell == "K562" & loop_data$factor == "RAD21"),c(1:4,6:7)], cell=rep("K562", length(which(loop_data$cell == "K562" & loop_data$factor == "RAD21")))))
annotated_loop_data.20k = rbind(annotated_loop_data.20k, cbind(loop_data[which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C"),c(1:4,6:7)], cell=rep("CH12", length(which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C")))))
annotated_loop_data.20k = add_te_intersections(annotated_loop_data.20k, loop_te_data)
annotated_loop_data.20k = add_motif_orientations(annotated_loop_data.20k, loop_motif_data)
annotated_loop_data.20k = intersect_cons_te_2(loop_conservation_data.20k, annotated_loop_data.20k)
dat_upset =  do_upset_plot(annotated_loop_data.s25k, "CH12", "GM12878", "20k")
pdf("upset_CH12-to-GM12878.20k.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()
dat.sets =  do_bar_plot(annotated_loop_data.s25k, "CH12", "GM12878", "20k")
pdf("stacked-bar_CH12-to-GM12878.20k.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

# 50k
loop_conservation_data.50k = read.table("../loop_conservation/RAD21_loops.all.50k.dat", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names=c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "IAB", "FDR", "mapped_chr_l", "mapped_start_l", "mapped_end_l", "qstrand_l", "mapped_chr_r", "mapped_start_r", "mapped_end_r", "qstrand_r", "id_l", "anchor_l", "chrom_l", "start_l", "end_l",  "id_r", "anchor_r", "chrom_r", "start_r", "end_r", "class", "cell_q", "cell_t"))
annotated_loop_data.50k = cbind(loop_data[which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"),c(1:4,6:7)], cell=rep("GM12878", length(which(loop_data$cell == "GM12878" & loop_data$factor == "RAD21"))))
annotated_loop_data.50k = rbind(annotated_loop_data.50k, cbind(loop_data[which(loop_data$cell == "K562" & loop_data$factor == "RAD21"),c(1:4,6:7)], cell=rep("K562", length(which(loop_data$cell == "K562" & loop_data$factor == "RAD21")))))
annotated_loop_data.50k = rbind(annotated_loop_data.50k, cbind(loop_data[which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C"),c(1:4,6:7)], cell=rep("CH12", length(which(loop_data$cell == "CH12" & loop_data$factor == "Hi-C")))))
annotated_loop_data.50k = add_te_intersections(annotated_loop_data.50k, loop_te_data)
annotated_loop_data.50k = add_motif_orientations(annotated_loop_data.50k, loop_motif_data)
annotated_loop_data.50k = intersect_cons_te_2(loop_conservation_data.50k, annotated_loop_data.50k)
dat_upset =  do_upset_plot(annotated_loop_data.s25k, "CH12", "GM12878", "50k")
pdf("upset_CH12-to-GM12878.50k.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()
dat.sets =  do_bar_plot(annotated_loop_data.s25k, "CH12", "GM12878", "50k")
pdf("stacked-bar_CH12-to-GM12878.50k.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()


## Figure 6D-F shows the contributions of enriched TE families
# to each conservation class for all species comparisons.

# Mouse to Human data (6D)
tmp = dat.CH12[which(dat.CH12$te_derived == TRUE),]
test = compile_class_enrichment_summary_data(tmp, "GM12878", "K562", useNA=TRUE)
test[is.na(test)] = 0
test$Human_Shared = test[,2]+test[,4]
test = test[,c(1,3,5,6)]
dat = melt(test)
dat$variable = factor(dat$variable, levels=c("Mouse", "Human_Shared", "No_Enr"))

pdf("TE-enrichment-class-by-conservation-class.mouse-to-human.3.pdf")
ggplot(dat, aes(x=class, y=value, fill=variable)) +
geom_bar(stat="identity", position="fill")
dev.off()
pdf("TE-enrichment-class-by-conservation-class.dotplot.mouse-to-human.3.pdf")
ggplot(dat, aes(x=class, y=value)) +
geom_point()
dev.off()

# Human to Mouse (6E)
tmp = dat.GM12878[which(dat.GM12878$te_derived == TRUE),]
test = compile_class_enrichment_summary_data(tmp, "CH12", useNA=TRUE)
tmp = dat.K562[which(dat.K562$te_derived == TRUE),]
test1 = compile_class_enrichment_summary_data(tmp, "CH12", useNA=TRUE)
test[,2:5] = test[,2:5] + test1[,2:5]
test = test[,c(1,2,4,5)]
dat = melt(test)
dat$variable = factor(dat$variable, levels=c("Human", "Shared", "No_Enr"))

pdf("TE-enrichment-class-by-conservation-class.human-to-mouse.pdf")
ggplot(dat, aes(x=class, y=value, fill=variable)) +
geom_bar(stat="identity", position="fill")
dev.off()
pdf("TE-enrichment-class-by-conservation-class.dotplot.human-to-mouse.pdf")
ggplot(dat, aes(x=class, y=value)) +
geom_point()
dev.off()

# Human to Human (6F)
tmp = dat.GM12878[which(dat.GM12878$te_derived == TRUE),]
test = compile_class_enrichment_summary_data(tmp, "K562", useNA=TRUE)
tmp = dat.K562[which(dat.K562$te_derived == TRUE),]
test1 = compile_class_enrichment_summary_data(tmp, "GM12878", useNA=TRUE)
test[,2:5] = test[,2:5] + test1[,2:5]
test = test[,c(1,2,4,5)]
dat = melt(test)
dat$variable = factor(dat$variable, levels=c("Human", "Shared", "No_Enr"))

pdf("TE-enrichment-class-by-conservation-class.human-to-human.2.pdf")
ggplot(dat, aes(x=class, y=value, fill=variable)) +
geom_bar(stat="identity", position="fill")
dev.off()
pdf("TE-enrichment-class-by-conservation-class.dotplot.human-to-human.2.pdf")
ggplot(dat, aes(x=class, y=value)) +
geom_point()
dev.off()


## Figure 6G compares TE age distributions across conservation classes for mouse-human and human-mouse comparisons.

# Mouse-human
tmp = dat.CH12
# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"GM12878_l_maps"] == FALSE) | (tmp$te_right & tmp[,"GM12878_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp
# Combine cell references
tmp1 = tmp[which(tmp$te_derived == TRUE),c("class_GM12878", "estAge")]
colnames(tmp1) = c("class", "estAge")
tmp2 = tmp[which(tmp$te_derived == TRUE),c("class_K562", "estAge")]
colnames(tmp2) = c("class", "estAge")
tmp1 = rbind(tmp1, tmp2)
# Plot the data
pdf("te-ages-by-cons-class.mouse-to-human.pdf")
ggplot(tmp1, aes(x=class, y=estAge/1000000)) +
geom_boxplot(notch=T)
dev.off()

# Human-Human
pdf("te-ages-by-cons-class.human-to-human.pdf")
# class_t2 always corresponds to CH12 target for human query cells
ggplot(annotated_loop_data[which(annotated_loop_data$te_derived == TRUE & (annotated_loop_data$cell_q == "GM12878" | annotated_loop_data$cell_q == "K562")),], aes(x=class_t2, y=estAge/1000000)) +
geom_boxplot(notch=T)
dev.off()

# Human-Mouse
pdf("te-ages-by-cons-class.human-to-mouse.pdf")
# class_t1 always corresponds to one human target cell or the other
ggplot(annotated_loop_data[which(annotated_loop_data$te_derived == TRUE & (annotated_loop_data$cell_q == "GM12878" | annotated_loop_data$cell_q == "K562")),], aes(x=class_t1, y=estAge/1000000)) +
geom_boxplot(notch=T)
dev.off()
