annotate_anchor_tss = function(dat, genesets=c("hg19"="/data/projects/adadiehl/genes/knownGenes.hg19.tss.sorted.bed","mm9"="/data/projects/adadiehl/genes/knownGenes.mm9.tss.sorted.bed"), cells=c("GM12878", "K562", "CH12"), species_index=c("CH12"="mm9", "GM12878"="hg19", "K562"="hg19")) {
    # Annotate loop anchors with their nearest TSS.
    ret = dat
    ret[,"tss_dist.r"] = ret[,"nearest_gene.r"] = ret[,"tss_dist.l"] = ret[,"nearest_gene.l"] = NA
    for (cell in cells) {
	genes = genesets[species_index[cell]]
        recs = annotate_tss(dat[which(dat$cell_q == cell), c("chrom1", "start1", "end1", "id_cp")], genes)
	ret[which(dat$cell_q == cell),"nearest_gene.l"] = recs[,"nearest_gene"]
	ret[which(dat$cell_q == cell),"tss_dist.l"] = recs[,"tss_dist"]
	recs = annotate_tss(dat[which(dat$cell_q == cell), c("chrom1", "start2", "end2", "id_cp")], genes)
        ret[which(dat$cell_q == cell),"nearest_gene.r"] = recs[,"nearest_gene"]
        ret[which(dat$cell_q == cell),"tss_dist.r"] = recs[,"tss_dist"]
    }
    return(ret)
}

annotate_tss = function (dat, genes) {
    # Annotate a set of BED coordinates with their nearest TSS. (based on the midpoint of the BED region)
    ret = dat
    ret$tss_dist = ret$nearest_gene = NA
    tmp_f = paste(paste(sample(0:9, 20, replace=TRUE), collapse=""), ".tmp", sep="")
    midpoints = dat[,2] + round( (dat[,3] - dat[,2]) / 2)
    ints = data.frame(chrom=dat[,1], start=midpoints, end=midpoints+1, name=dat[,4], stringsAsFactors=FALSE)
    ints = ints[order(ints[,1], ints[,2]),]
    write.table(ints, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", file=tmp_f)
    cmd_args = c("closest", "-a", tmp_f, "-b", genes)
    res = system2("bedtools", args=cmd_args, stdout=TRUE, stderr="")
    unlink(tmp_f)
    for (i in 1:length(res)) {
        rec = unlist(strsplit(res[i],"\t"))
	ret[which(ret[,4] == rec[4]), "nearest_gene"] = rec[8]
	ret[which(ret[,4] == rec[4]), "tss_dist"] = as.numeric(rec[2]) - as.numeric(rec[6])
    }
    return(ret)
}

add_expr_fc = function(dat, expr, species_index=c("CH12"="mm9", "GM12878"="hg19", "K562"="hg19")) {
    # Add expression fold-changes for nearest genes to each anchor in comparison to target cell types.
    ret = dat
    ret$expr_t2.r = ret$expr_t2.l = ret$expr_t1.r = ret$expr_t1.l = ret$expr_q.r = ret$expr_q.l = NA
    ret$expr_fc_t2.r = ret$expr_fc_t2.l = ret$expr_fc_t1.r = ret$expr_fc_t1.l = NA
    for (i in 1:nrow(ret)) {
        t1 = ret[i,"cell_t1"]
	t2 = ret[i,"cell_t2"]
	gene_l = ret[i,"nearest_gene.l"]
	gene_r = ret[i,"nearest_gene.r"]
	gene_name_col =	paste("gene", species_index[ret[i,"cell_q"]], sep="_")
	print(gene_name_col)
	print(gene_l)
	print(gene_r)

	# Nearest gene to left anchor
	e = expr[which(expr[,gene_name_col] == gene_l),]
	if (nrow(e) > 0) {
	    print(e)
	    ret[i,"expr_q.l"] = e[ret[i,"cell_q"]]
	    ret[i,"expr_t1.l"] = e[ret[i,"cell_t1"]]
	    ret[i,"expr_t2.l"] = e[ret[i,"cell_t2"]]
	    ret[i,"expr_fc_t1.l"] = e[ret[i,"cell_q"]] / e[ret[i,"cell_t1"]]
	    ret[i,"expr_fc_t2.l"] = e[ret[i,"cell_q"]] / e[ret[i,"cell_t2"]]
	}

	# Nearest gene to right anchor
	e = expr[which(expr[,gene_name_col] == gene_r),]
	if (nrow(e) > 0) {
            ret[i,"expr_q.r"] = e[ret[i,"cell_q"]]
            ret[i,"expr_t1.r"] = e[ret[i,"cell_t1"]]
            ret[i,"expr_t2.r"] = e[ret[i,"cell_t2"]]
            ret[i,"expr_fc_t1.r"] = e[ret[i,"cell_q"]] / e[ret[i,"cell_t1"]]
            ret[i,"expr_fc_t2.r"] = e[ret[i,"cell_q"]] / e[ret[i,"cell_t2"]]
	}
    }
    return(ret)
}

