# Assembles data from loop cross-mapping into data tables.

## Load in commands
source("intersect_cons_te.R")

## Read in data

# Loop conservation classification results
loop_conservation_data = read.table("../loop_conservation/RAD21_loops.all.dat", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(loop_conservation_data) = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "IAB", "FDR", "mapped_chr_l", "mapped_start_l", "mapped_end_l", "qstrand_l", "mapped_chr_r", "mapped_start_r", "mapped_end_r", "qstrand_r", "id_l", "anchor_l", "chrom_l", "start_l", "end_l",  "id_r", "anchor_r", "chrom_r", "start_r", "end_r", "class", "cell_q", "cell_t")

# TE-Loop intersections
dat.CH12 = read.table("../te_loop_overlaps/data.CH12.txt", stringsAsFactors=FALSE)
dat.K562 = read.table("../te_loop_overlaps/data.K562.txt", stringsAsFactors=FALSE)
dat.GM12878 = read.table("../te_loop_overlaps/data.GM12878.txt", stringsAsFactors=FALSE)


## Intersect the conservation data with the TE-loop intersections
dat.GM12878 = intersect_cons_te(loop_conservation_data, dat.GM12878, "CH12", "K562")
dat.K562 = intersect_cons_te(loop_conservation_data, dat.K562, "CH12", "GM12878")
dat.CH12 = intersect_cons_te(loop_conservation_data, dat.CH12, "GM12878", "K562")


## Compile tabular data to incorporate into supplementary table 7
# First put all cell-wise data into a single list...  
dat = mget(ls(pattern = "dat.[A-Z]"))
names(dat) = lapply(names(dat), function (e) {unlist(strsplit(e, '.', fixed=TRUE))[2]}) 
write.table(compile_summary_data(dat), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, file="TE-derived-loop-cons.txt")


## Generate stacked bar charts for supplementary figure 5A
dat.classes = read.table("RAD21_loops.stats.dat", header=TRUE, stringsAsFactors=FALSE, sep = "\t")
dat.classes$B2 = dat.classes$B2S+dat.classes$B2D
dat.classes = dat.classes[,c(1,5,10,7,6,3,4,2)]

library(ggplot2)
library(reshape)
pdf("class-breakdown-by-cell-comparison.scaled.pdf")
ggplot(melt(dat.classes), aes(x=File, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()


## Generate the species-comparison-wise stacked bar charts for Figure 3B
tmp = data.frame(matrix(unlist(c((dat.classes[1,2:8] + dat.classes[2,2:8]), (dat.classes[3,2:8] + dat.classes[5,2:8]), (dat.classes[4,2:8] + dat.classes[6,2:8])), nrow=3, byrow=T))
colnames(tmp) = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")
tmp$Comparison =  c("mouse-human", "human-mouse", "human-human")

pdf("class-breakdown-by-species-comparison.scaled.pdf")
ggplot(melt(tmp), aes(x=Comparison, y=value, fill=variable)) +
geom_bar(stat='identity', position='fill')
dev.off()


## Figures 3C-3D, and supplementary figure 7 are a combination of bar plots and Upset plots. The bar chart at the top of the figure is plotted separately from the set size bar chart and set illustration dot plot, and then combined by hand in illustrator to produce the final figure...

# Generate the bar plot components...

# CH12 to GM12878
dat.sets = list("Conserved" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B2" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B2"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B2" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B1" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B1"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B1" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B0" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B0"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "B0" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "N1A" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1A"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1A" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1B"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N1B" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0"),]), nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),]))
               )
dat.sets = as.data.frame(t(as.data.frame(dat.sets)))
dat.sets[,1] = dat.sets[,1] - dat.sets[,2]
colnames(dat.sets) = c("No TE", "TE")
dat.sets$cat = factor(rownames(dat.sets), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("stacked-bar_CH12-to-GM12878.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

 # CH12 to K562
dat.sets = list("Conserved" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "C"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "C" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B2" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "B2"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "B2" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B1" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "B1"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "B1" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "B0" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "B0"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "B0" & (dat.CH12$te_left | dat.CH12$te_right)),])),
                "N1A" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "N1A"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "N1A" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "N1B"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "N1B" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.CH12[which(dat.CH12$class_K562 == "N0"),]), nrow(dat.CH12[which(dat.CH12$class_K562 == "N0" & ( (dat.CH12$te_left & dat.CH12$K562_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$K562_r_maps == FALSE) )),]))
               )
dat.sets = as.data.frame(t(as.data.frame(dat.sets)))
dat.sets[,1] = dat.sets[,1] - dat.sets[,2]
colnames(dat.sets) = c("No TE", "TE")
dat.sets$cat = factor(rownames(dat.sets), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("stacked-bar_CH12-to-K562.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()

# GM12878 to K562
dat.sets = list("Conserved" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "C"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "C" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B2" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B2"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B2" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B1" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B1"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B1" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "B0" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B0"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "B0" & (dat.GM12878$te_left | dat.GM12878$te_right)),])),
                "N1A" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1A"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1A" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1B"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N1B" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N0"),]), nrow(dat.GM12878[which(dat.GM12878$class_K562 == "N0" & ( (dat.GM12878$te_left & dat.GM12878$K562_l_maps == FALSE) | (dat.GM12878$te_right & dat.GM12878$K562_r_maps == FALSE) )),]))
               )
dat.sets = as.data.frame(t(as.data.frame(dat.sets)))
dat.sets[,1] = dat.sets[,1] - dat.sets[,2]
colnames(dat.sets) = c("No TE", "TE")
dat.sets$cat = factor(rownames(dat.sets), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("stacked-bar_GM12878-to-K562.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()
 
# K562 to GM12878
dat.sets = list("Conserved" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "C"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "C" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B2" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "B2"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "B2" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B1" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "B1"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "B1" & (dat.K562$te_left | dat.K562$te_right)),])),
                "B0" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "B0"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "B0" & (dat.K562$te_left | dat.K562$te_right)),])),
                "N1A" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "N1A"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "N1A" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),])),
                "N1B" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "N1B"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "N1B" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),])),
                "N0" = c(nrow(dat.K562[which(dat.K562$class_GM12878 == "N0"),]), nrow(dat.K562[which(dat.K562$class_GM12878 == "N0" & ( (dat.K562$te_left & dat.K562$GM12878_l_maps == FALSE) | (dat.K562$te_right & dat.K562$GM12878_r_maps == FALSE) )),]))
               )
dat.sets = as.data.frame(t(as.data.frame(dat.sets)))
dat.sets[,1] = dat.sets[,1] - dat.sets[,2]
colnames(dat.sets) = c("No TE", "TE")
dat.sets$cat = factor(rownames(dat.sets), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))
 
pdf("stacked-bar_K562-to-GM12878.pdf")
ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
geom_bar(stat='identity', position='fill')
dev.off()


# Generate the upset plot components...
library(UpSetR)

# CH12 to GM12878
dat_upset = c("Conserved" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "C"),]),
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
              "N0&TE" = nrow(dat.CH12[which(dat.CH12$class_GM12878 == "N0" & ( (dat.CH12$te_left & dat.CH12$GM12878_l_maps == FALSE) | (dat.CH12$te_right & dat.CH12$GM12878_r_maps == FALSE) )),])

pdf("upset_CH12-to-GM12878.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()

# CH12 to K562
dat_upset = c("Conserved" = nrow(dat.CH12[which(dat.CH12$class_K562 == "C"),]),
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

pdf("upset_CH12-to-K562.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()


# K562 to GM12878
dat_upset = c("Conserved" = nrow(dat.K562[which(dat.K562$class_GM12878 == "C"),]),
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

pdf("upset_K562-to-GM12878.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()

# GM12878 to K562
dat_upset = c("Conserved" = nrow(dat.GM12878[which(dat.GM12878$class_GM12878 == "C"),]),
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

pdf("upset_GM12878-to-K562.pdf")
upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
dev.off()


## Supplementary figure 5B-G: Loop strength relative to loop conservation

# Box plot components:

# GM21878 query
tmp = dat.GM12878
tmp$IAB = loop_conservation_data[which(loop_conservation_data$id == dat.GM12878$id_cp & loop_conservation_data$cell_t == "CH12"),"IAB"]
# Does not matter which cell_t you use here
tmp$class_CH12 = factor(tmp$class_CH12, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
tmp$class_K562 = factor(tmp$class_K562, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

pdf("IAB-boxplots.GM12878-to-CH12.pdf")
ggplot(tmp[which(  tmp$IAB <= quantile(tmp$IAB, probs=c(0.9)) ),], aes(x=class_CH12, y=IAB, fill=te_derived)) +
geom_boxplot(notch=T)
dev.off()

pdf("IAB-boxplots.GM12878-to-K562.pdf")
gplot(tmp[which(  tmp$IAB <= quantile(tmp$IAB, probs=c(0.9)) ),], aes(x=class_K562, y=IAB, fill=te_derived)) +
geom_boxplot(notch=T)
dev.off()

# CH12 Query
tmp = dat.CH12
tmp$IAB = loop_conservation_data[which(loop_conservation_data$id %in% dat.CH12$id_cp & loop_conservation_data$cell_t == "GM12878"),"IAB"]
tmp$class_GM12878 = factor(tmp$class_GM12878, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
tmp$class_K562 = factor(tmp$class_K562, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"GM12878_l_maps"] == FALSE) | (tmp$te_right & tmp[,"GM12878_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

pdf("IAB-boxplots.CH12-to-GM12878.pdf")
ggplot(tmp[which(  tmp$IAB <= quantile(tmp$IAB, probs=c(0.9)) ),], aes(x=class_GM12878, y=IAB, fill=te_derived)) +
geom_boxplot(notch=T)
dev.off()

pdf("IAB-boxplots.CH12-to-K562.pdf")
ggplot(tmp[which(  tmp$IAB <= quantile(tmp$IAB, probs=c(0.9)) ),], aes(x=class_K562, y=IAB, fill=te_derived)) +
geom_boxplot(notch=T)
dev.off()

# K562 Query
tmp = dat.K562
tmp$IAB = loop_conservation_data[which(loop_conservation_data$id %in% dat.K562$id_cp & loop_conservation_data$cell_t == "CH12"),"IAB"]
tmp$class_GM12878 = factor(tmp$class_GM12878, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
tmp$class_CH12 = factor(tmp$class_CH12, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

pdf("IAB-boxplots.K562-to-GM12878.pdf")
ggplot(tmp[which(  tmp$IAB <= quantile(tmp$IAB, probs=c(0.9)) ),], aes(x=class_GM12878, y=IAB, fill=te_derived)) +
geom_boxplot(notch=T)
dev.off()

pdf("IAB-boxplots.K562-to-CH12.pdf")
ggplot(tmp[which(  tmp$IAB <= quantile(tmp$IAB, probs=c(0.9)) ),], aes(x=class_CH12, y=IAB, fill=te_derived)) +
geom_boxplot(notch=T)
dev.off()

# heat map components:

source("loop_score_wilcox_tests.R")

# GM12878 to CH12
tmp = dat.GM12878
tmp$IAB = loop_conservation_data[which(loop_conservation_data$id == dat.GM12878$id_cp & loop_conservation_data$cell_t == "CH12"),"IAB"]
# Does not matter which cell_t you use here
tmp$class_CH12 = factor(tmp$class_CH12, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
tmp$class_K562 = factor(tmp$class_K562, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

dat = build_wilcox_matrix(tmp, "CH12")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-score-comparisons_p-vals-heat.GM12878-to-CH12.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()

# GM12878 to K562
dat = build_wilcox_matrix(tmp, "K562")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-score-comparisons_p-vals-heat.GM12878-to-K562.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()

# CH12 to GM12878
tmp = dat.CH12
tmp$IAB = loop_conservation_data[which(loop_conservation_data$id %in% dat.CH12$id_cp & loop_conservation_data$cell_t == "GM12878"),"IAB"]
tmp$class_GM12878 = factor(tmp$class_GM12878, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
tmp$class_K562 = factor(tmp$class_K562, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"GM12878_l_maps"] == FALSE) | (tmp$te_right & tmp[,"GM12878_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

dat = build_wilcox_matrix(tmp, "GM12878")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-score-comparisons_p-vals-heat.CH12-to-GM12878.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()

# CH12 to K562
dat = build_wilcox_matrix(tmp, "K562")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-score-comparisons_p-vals-heat.CH12-to-K562.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()

# K562 to GM12878
tmp = dat.K562
tmp$IAB = loop_conservation_data[which(loop_conservation_data$id %in% dat.K562$id_cp & loop_conservation_data$cell_t == "CH12"),"IAB"]
tmp$class_GM12878 = factor(tmp$class_GM12878, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
tmp$class_CH12 = factor(tmp$class_CH12, levels = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

dat = build_wilcox_matrix(tmp, "GM12878")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-score-comparisons_p-vals-heat.K562-to-GM12878.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()

# K562 to CH12
dat = build_wilcox_matrix(tmp, "CH12")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-score-comparisons_p-vals-heat.K562-to-CH12.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()


## Supplementary figure 6: comparison of TE-derived and non-TE CTCF ChIP-seq score distributions in each cell
all_ctcf = read.table("all_ctcf_peaks.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)

pdf("ctcf-signal_TE-v-nonTE.all.pdf")
ggplot(all_ctcf[which(all_ctcf$signalValue <= quantile(all_ctcf$signalValue, probs=c(0.9))),], aes(x=signalValue, fill=is_te))+
geom_density(alpha=0.3)
dev.off()

pdf("ctcf-signal_TE-v-nonTE.GM12878.pdf")
ggplot(all_ctcf[which(all_ctcf$cell == "GM12878"),], aes(x=signalValue, fill=is_te))+
geom_density(alpha=0.3)
dev.off()

pdf("ctcf-signal_TE-v-nonTE.K562.pdf")
ggplot(all_ctcf[which(all_ctcf$cell == "K562"),], aes(x=signalValue, fill=is_te))+
geom_density(alpha=0.3)
dev.off()

pdf("ctcf-signal_TE-v-nonTE.CH12.pdf")
ggplot(all_ctcf[which(all_ctcf$cell == "CH12" & all_ctcf$signalValue <= 100),], aes(x=signalValue, fill=is_te))+
geom_density(alpha=0.3)
dev.off()


## Figure 3G compares TE age distributions across conservation classes for mouse-human and human-mouse comparisons.

# Mouse-human comparison is straightforward...
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

# Generate the heat map of p-values for all pairwise Wilcoxon rank-sum tests.
dat = wilcox_age_matrix(dat.CH12, "GM12878", "K562", alt="g")
dat.m = melt(dat)
dat.m$X1 = factor(dat.m$X1, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
dat.m$X2 = factor(dat.m$X2, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))

pdf("class-age-comparisons_p-vals-heat.mouse-to-human.pdf")
ggplot(dat.m, aes(x=X2, y=X1)) +
geom_tile(aes(fill=value), colour="black") +
coord_fixed(ratio=1) +
scale_fill_gradient(high="white", low="#6F2E91", limits=c(0,1))
dev.off()

# Test for linear correlation between conservation classes and TE age in the mouse-human comparison
tmp1$class = as.numeric(tmp1$class)

# Scatter plot with lm line overlayed.
ggplot(tmp1, aes(x=class, y=estAge)) +
geom_point() +
geom_smooth(method="lm")

# Correlation stats...
summary(lm(class ~ estAge, data=tmp1))


# Human-human comparison...
tmp = dat.GM12878

# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp

# Combine cell references
tmp2 = tmp[which(tmp$te_derived == TRUE),c("class_K562", "estAge")]
colnames(tmp2) = c("class", "estAge")

# Get the K562 data
tmp = dat.K562

# Separate out TE-derived loops
tmp1 = tmp[which( ( ( (tmp$te_left & tmp[,"CH12_l_maps"] == FALSE) | (tmp$te_right & tmp[,"CH12_r_maps"] == FALSE) ) ) ) ,]
tmp1 = rbind(tmp1, tmp[which( !(tmp$id_cp %in% tmp1$id_cp) & (tmp$te_left | tmp$te_right) ),])
tmp$te_derived = tmp$id_cp %in% tmp1$id_cp
tmp1 = tmp[which(tmp$te_derived == TRUE),c("class_GM12878", "estAge")]
colnames(tmp1) = c("class", "estAge")

# Roll in the data
tmp1 = rbind(tmp1, tmp2)

# Plot the data
pdf("te-ages-by-cons-class.human-to-human.pdf")
ggplot(tmp1, aes(x=class, y=estAge/1000000)) +
geom_boxplot(notch=T)
dev.off()

# No Wilcoxon p-value matrix here. Do the test for linear correlation, though, for comparison.
tmp1$class = as.numeric(tmp1$class)

ggplot(tmp1, aes(x=class, y=estAge)) +
geom_point() +
geom_smooth(method="lm")

summary(lm(class ~ estAge, data=tmp1))


## Figure 3E shows the contributions of enriched TE families to each conservation class in mouse-human, and 3F in  human-mouse comparisons.
# Mouse to Human data (3E)
tmp = dat.CH12[which(dat.CH12$te_name_l %in% enr_te_fams$name | dat.CH12$te_name_r %in% enr_te_fams$name & dat.CH12$te_derived == TRUE),]
test = compile_class_enrichment_summary_data(tmp, "GM12878", "K562")
# Combine shared and human-enriched categories for clarity
test[is.na(test)] = 0
test$Human_Shared = test[,2]+test[,4]
test = test[,c(1,3,5)]

pdf("TE-enrichment-class-by-conservation-class.mouse-to-human.pdf")
ggplot(melt(test), aes(x=class, y=value, fill=variable)) +
geom_bar(stat="identity", position="fill")
dev.off()

# Dot plot to indicate set sizes -- this will be hand-overlayed onto the bar plot.
pdf("TE-enrichment-class-by-conservation-class.dotplot.mouse-to-human.pdf")
ggplot(melt(test), aes(x=class, y=value)) +
geom_point()
dev.off()

# Human to Human (3F)
tmp = dat.GM12878[which(dat.GM12878$te_name_l %in% enr_te_fams$name | dat.GM12878$te_name_r %in% enr_te_fams$name & dat.GM12878$te_derived == TRUE),]
test = compile_class_enrichment_summary_data(tmp, "K562")
tmp = dat.K562[which(dat.K562$te_name_l %in% enr_te_fams$name | dat.K562$te_name_r %in% enr_te_fams$name & dat.K562$te_derived == TRUE),]
test1 = compile_class_enrichment_summary_data(tmp, "GM12878")
test[,2:4] = test[,2:4] + test1[,2:4]

pdf("TE-enrichment-class-by-conservation-class.human-to-human.pdf")
ggplot(melt(test), aes(x=class, y=value, fill=variable)) +
geom_bar(stat="identity", position="fill")
dev.off()

# Dot plot to indicate set sizes -- this will be hand-overlayed onto the bar plot.
pdf("TE-enrichment-class-by-conservation-class.dotplot.human-to-human.pdf")
ggplot(melt(test), aes(x=class, y=value)) +
geom_point()
dev.off()


## Figure 3H shows the spatial distribution of phastCons conservation scores across 500 bp
## windows surrounding the annotated CTCF ChIP-seq peak in TE-derived loop anchors, broken
## down by conservation class.

# Human-mouse data are used in Figure 3H:
test2 = compile_class_bigwig_scores(dat.GM12878, loop_te_data, "/data/UCSC/PHASTCONS/hg19/phastCons46way.placental.bigWig", "CH12", size=1000)
# Convert data to long format matrix. For some reason, the "melt" function doesn't play well here...
dat2 = convert_to_long(test2)
# Plot the data as lines...
pdf("loop-conservation_vs_phastcons.GM12878-to-CH12.pdf")
ggplot(dat2, aes(x=position, y=value, colour=cat))+
geom_line()
dev.off()

# Mouse-human shows the same trends, but much noisier due to the small set sizes.
test = compile_class_bigwig_scores(dat.CH12, hic_te_data, "/data/UCSC/PHASTCONS/mm9/phastCons30way.placental.bigWig", "GM12878", "K562", size=1000)
dat = convert_to_long(test)
ggplot(dat, aes(x=position, y=value, colour=cat))+
geom_line()

# Human-human data have similar distributions for all classes, with only slight conservation drop at the peaks.
test3 = compile_class_bigwig_scores(dat.GM12878, loop_te_data, "/data/UCSC/PHASTCONS/hg19/phastCons46way.placental.bigWig", "K562", size=1000, cats=c("C","B2","B1","B0"))
dat3 = convert_to_long(test3, cats=c("C","B2","B1","B0"))
ggplot(dat3, aes(x=position, y=value, colour=cat))+
geom_line()
