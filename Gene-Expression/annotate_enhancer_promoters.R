add_histmod = function(dat, bigwig_path, mod, chrom_col, start_col, end_col, cells=c("GM12878", "K562", "CH12"), stat="mean", na.rm=TRUE) {
    # Add annotations for a single histone modification across all cells.
    ret = data.frame(cell=dat$cell_q, res=rep(NA, nrow(dat)), stringsAsFactors=FALSE, row.names=rownames(dat))
    for (cell in cells) {
    	bigwig = paste(bigwig_path, paste(cell, mod, "bigWig", sep="."), sep="/")
        ret[which(dat$cell_q == cell),"res"] = get_bigwig_scores(bed_frame(dat[which(dat$cell_q == cell),], chrom_col, start_col, end_col), bigwig, stat=stat, na.rm=na.rm)	
    }
    return(ret$res)
}

add_anchor_histmods = function(dat, bigwig_path, mod, cells=c("GM12878", "K562", "CH12"), stat="mean", na.rm=TRUE) {
    # Add histone mods for right and left anchors
    ret = dat
    ret[,paste(mod, stat, "l", sep=".")] = add_histmod(dat, bigwig_path, mod, "chrom1", "start1", "end1", cells=cells, stat=stat, na.rm=na.rm)
    ret[,paste(mod, stat, "r", sep=".")] = add_histmod(dat, bigwig_path, mod, "chrom1", "start2", "end2", cells=cells, stat=stat, na.rm=na.rm)
    return(ret)
}

add_all_histmods = function(dat, mods=c("H3K4me3", "H3K4me1", "H3K9me3", "H3K27me3", "H3K27ac"), bigwig_path="../gene_expression", cells=c("GM12878", "K562", "CH12"), stat="max", na.rm=TRUE) {
    # Add histone mod annotations for the given set of mods.
    ret = dat
    for (mod in mods) {
        ret = add_anchor_histmods(ret, bigwig_path, mod, cells=cells, stat=stat, na.rm=na.rm)
    }
    return(ret)
}

rescale = function(dat) {
    # Perform a simple min-max scaling to a range of values.
    ret = unlist(lapply(dat, function(x) { return( (x-min(dat, na.rm=TRUE)) / (max(dat, na.rm=TRUE)-min(dat, na.rm=TRUE)) ) }))
}

standardize = function(dat) {
    # Perform a z-score standardization to a range of values.
    ret = unlist(lapply(dat, function(x) { return( (x - mean(dat, na.rm=TRUE)) / sd(dat, na.rm=TRUE) ) }))
    return(ret)
}

norm_all_histmods = function(dat, mods=c("H3K4me3", "H3K4me1", "H3K9me3", "H3K27me3", "H3K27ac"), cells=c("GM12878", "K562", "CH12"), norm="rescale", stat="max") {
    # Normalize values using the given method for all samples.
    ret = dat
    for (mod in mods) {
        ret = norm_anchor_histmods(ret, mod, cells=cells, norm=norm, stat=stat)
    }
    return(ret)
}

norm_anchor_histmods = function(dat, mod, cells=c("GM12878", "K562", "CH12"), norm="rescale", stat="max") {
    # Normalize histone mods for right and left anchors
    ret = dat
    ret[,paste(mod, stat, norm, "l", sep=".")] = norm_histmod(dat, mod, col=paste(mod, stat, "l", sep="."), cells=cells, norm=norm)
    ret[,paste(mod, stat, norm, "r", sep=".")] = norm_histmod(dat, mod, col=paste(mod, stat, "r", sep="."), cells=cells, norm=norm)
    return(ret)
}

norm_histmod = function(dat, mod, col, cells=c("GM12878", "K562", "CH12"), norm="rescale") {
    # Normalize annotations for a single histone modification across all cells.
    ret = data.frame(cell=dat$cell_q, res=rep(NA, nrow(dat)), stringsAsFactors=FALSE, row.names=rownames(dat))
    for (cell in cells) {
        ret[which(dat$cell_q == cell),"res"] = match.fun(norm)(dat[which(dat$cell_q == cell),col])
    }
    return(ret$res)
}
