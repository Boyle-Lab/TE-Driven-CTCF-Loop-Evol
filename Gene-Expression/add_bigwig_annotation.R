get_bigwig_scores = function(dat, bigwig, stat="mean", na.rm=TRUE) {
    # Uses rtracklayer to retrieve values from a bigwig
    # file using the given stat to aggregate over all
    # columns. Valid stats are "mean", "median", "max",
    # "min", "stdev", "sum", and "full".
    # "dat" = Data frame including "chrom", "start", and "end"
    #         columns defining genomic regions to intersect
    #	      with bigWig scores.
    # "bigwig" = Path to a bigWig file.
    # "stat" = the statistic to return for each region.

    require("rtracklayer")
    # Get full data from the bigWig file
    res = import.bw(con = BigWigFile(bigwig),
		    selection = makeGRangesFromDataFrame(dat,
							 ignore.strand=TRUE),
		    as = "NumericList")
		    
    res = lapply(res, drop)

    if (stat == "mean") {
        return(unlist(lapply(res, function(x) { return(mean(x, na.rm=na.rm)) }) ))
    }
    if (stat == "sum") {
        return(unlist(lapply(res, function(x) { return(sum(x, na.rm=na.rm)) }) ))
    }
    if (stat == "median") {
        return(unlist(lapply(res, function(x) { return(median(x, na.rm=na.rm)) })))
    }
    if (stat == "max") {
        return(unlist(lapply(res, function(x) { return(max(x, na.rm=na.rm)) })))
    }
    if (stat == "min") {
	return(unlist(lapply(res, function(x) { return(min(x, na.rm=na.rm)) })))
    }
    if (stat == "stdev") {
	return(unlist(lapply(res, function(x) { return(sd(x, na.rm=na.rm)) })))
    }
    return(res)
}

bed_frame = function(dat, chrom_colname, start_colname, end_colname) {
    # Return a data frame with column name conventions for
    # granges bed: chrom, start, end; based on the given
    # colnames from dat.
    return(data.frame('chrom' = dat[,chrom_colname],
    	  	      'start' = dat[,start_colname],
		      'end' = dat[,end_colname],
		      row.names=rownames(dat),
		      stringsAsFactors = FALSE))
}
