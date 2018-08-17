intersect_cons_te = function(dat.cons, dat.te, cell1, cell2) {
    # Associates loop data/TE/motif intersection data with loop conservation class assignments.
    #
    # dat.te is a data frame containing TE/motif intersection data for a set of loops (typically
    # representing all annotated loops from a given cell type). Each row in this data frame is
    # a distinct loop, with rownames set to the loop ID number (This is used to perform the
    # intersection with conservation data.)
    #
    # dat.cons is a data frame containing conservation class assignments and associated loop
    # anchor overlap data for all loops in dat.te (and possibly others). All loops in dat.te
    # must be present in dat.cons, with a record for both cell types, or there will be errors.
    #
    # cell1 is the name of the first target cell type. This is the cell that will be used to
    # check for mappability, so if data are being mapped across species and only one cell is
    # from the target species, make sure to supply it as cell1!
    #
    # cell2 is the name of the second target cell type. Only class and overlap with the query
    # cell will be retrieved for this cell type, since mappability will be the same as with
    # the query cell versus cell1 and cell1 versus cell2. If cell1 and cell2 are from different
    # species, cell2 should be assigned to the cell from the same species as the query.
    
    ret = dat.te
    
    map_c1_l = paste(cell1, "l_maps", sep="_")
    map_c1_r = paste(cell1, "r_maps", sep="_")
    map_c2_l = paste(cell2, "l_maps", sep="_")
    map_c2_r = paste(cell2, "r_maps", sep="_")

    ret[,map_c2_r] = ret[,map_c2_l] = ret[,map_c1_r] = ret[,map_c1_l]= FALSE
    
    ol_c1_l = paste(cell1, "l_overlaps", sep="_")
    ol_c1_r = paste(cell1, "r_overlaps", sep="_")
    ol_c2_l = paste(cell2, "l_overlaps", sep="_")
    ol_c2_r = paste(cell2, "r_overlaps", sep="_")

    ret[,ol_c2_r] = ret[,ol_c2_l] = ret[,ol_c1_r] = ret[,ol_c1_l] = FALSE

    cl_c1 = paste("class", cell1, sep="_")
    row_c1 = paste("row", cell1, sep="_")
    cl_c2 = paste("class", cell2, sep="_")
    row_c2 = paste("row", cell2, sep="_")

    ret[,row_c2] = ret[,row_c1] = ret[,cl_c2] = ret[,cl_c1] = NA

    #print(head(ret))
    
    for (id in rownames(ret)) {
        x = dat.cons[which(dat.cons$id == id),]
	
	ret[which(ret$id == id), c(cl_c1, row_c1,
			           map_c1_l, map_c1_r,
				   ol_c1_l, ol_c1_r)] = c(x[which(x$cell_t == cell1),"class"],
			    	 		 	  rownames(x[which(x$cell_t == cell1),]),
			    	 		 	  (x[which(x$cell_t == cell1),"mapped_chr_l"] != "."),
							  (x[which(x$cell_t == cell1),"mapped_chr_r"] != "."),
							  (x[which(x$cell_t == cell1),"chrom_l"] != "."),
							  (x[which(x$cell_t == cell1),"chrom_r"] != "."))
	ret[which(ret$id == id), c(cl_c2, row_c2,
			           map_c2_l, map_c2_r,
		            	   ol_c2_l, ol_c2_r)] = c(x[which(x$cell_t == cell2),"class"],
			                                  rownames(x[which(x$cell_t == cell2),]),
						   	  (x[which(x$cell_t == cell2),"mapped_chr_l"] != "."),
                                                          (x[which(x$cell_t == cell2),"mapped_chr_r"] != "."),
							  (x[which(x$cell_t == cell2),"chrom_l"] != "."),
                                                          (x[which(x$cell_t == cell2),"chrom_r"] != "."))
    }
    return(ret)
}

compile_summary_data_inner = function(dat, cell) {
    # Counts the number of TE overlaps for non-mapping/non-overlapping anchors in each loopconservation class.
    
    ret = c("C" = 0, "B2" = 0, "B1" = 0, "B0" = 0, "N1A" = 0,"N1B" = 0, "N0" = 0)
    class_col = paste("class", cell, sep="_")
    map_col_l =	paste(cell, "l_maps", sep="_")
    map_col_r = paste(cell, "r_maps", sep="_")
    for (cat in dat[,class_col]) {
        if (cat == "C") {
	    ret[cat] = nrow(dat[which(dat[,class_col] == cat & (dat$te_left | dat$te_right)),])
	} else if (cat == "B0" || cat == "B1") {	
	    ret[cat] = nrow(dat[which(dat[,class_col] == cat & (dat$te_left | dat$te_right)),])
	} else if (cat == "B2") {
	    ret[cat] =  nrow(dat[which(dat[,class_col] == "B2" & (dat$te_left | dat$te_right)),])
	} else {
	    if (map_col_l %in% colnames(dat)) {
 	        ret[cat] = nrow(dat[which(dat[,class_col] == cat & ( (dat$te_left & dat[,map_col_l] == FALSE) | (dat$te_right & dat[,map_col_r] == FALSE) ) ),])
		
	    }
	}
    }
    return(ret)
}

compile_summary_data = function(dat, cells = c("GM12878", "CH12", "K562"), cats = c("C" = 0, "B2" = 0, "B1" = 0, "B0" = 0, "N1A" = 0,"N1B" = 0, "N0" = 0)) {
    # Compile a table of TE-conservation intersection data for all cell combinations.
    # "dat" is a list of data frames, with names set equal to the cells in the cells list.
    ret = as.data.frame(matrix(nrow = (length(cells)**2)-length(cells), ncol = length(cats)+2))
    colnames(ret) = c("cell1", "cell2", names(cats))
    i = 1
    for (cell1 in cells) {
        for (cell2 in cells) {
	    if (cell1 == cell2) {
	        next
	    }
	    print(c(cell1, cell2))
	    ret[i,"cell1"] = cell1
	    ret[i,"cell2"] = cell2
	    ret[i,3:ncol(ret)] = compile_summary_data_inner(dat[[ cell1 ]], cell2)
	    i = i+1
        }
    }
    return(ret)
}

add_te_ages = function(dat.loops, dat.te, rate) {
  # Add the percent divergence and estimated ages for overlapping TE
  dat = dat.loops
  dat$pctDiv = NA
  dat$estAge = NA
  for (id_cp in dat$id_cp) {          
        x = dat.te[which(dat.te$id_cp == id_cp),]
        if (nrow(x) > 0) {
            loop = dat[which(dat$id_cp == id_cp),]
	    #if (nrow(x) > 1) { print(x) }
	    # If there is more than 1 TE contributing to the loop, take the average.
	    pctDiv = 0
            for (row in 1:nrow(x)) {
                #print(x[row,])
		pctDiv = pctDiv + x[row,"pctdiv_rp"]
            }
	    pctDiv = pctDiv / nrow(x)
	    estAge = (pctDiv/100) / rate
	    dat[which(dat$id_cp == id_cp), "pctDiv"] = pctDiv
            dat[which(dat$id_cp == id_cp), "estAge"] = estAge
        }
    }
    return(dat)
}

compile_class_enrichment_summary_data = function(dat.loops, cell_1, cell_2=NULL, cats = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"), useNA=FALSE) {
    # Compile a data table summarizing the contribution of different
    # TE enrichment categories to the composition within each conservation class.
    class_col_1 = paste("class", cell_1, sep="_")
    if (!is.null(cell_2)) {
        class_col_2 = paste("class", cell_2, sep="_")
    }
    
    ret = as.data.frame(matrix(nrow=length(cats), ncol=4))
    colnames(ret) = c("class", "Human", "Mouse", "Shared")
    if (useNA) {
        ret$No_Enr = NA
    }
    for (i in 1:length(cats)) {
        cat = cats[i]
        ret[i,"class"] = cat

	tmp = dat.loops[which(dat.loops[,class_col_1] == cat),]
	if (!is.null(cell_2)) {	
	    tmp = rbind(tmp, dat.loops[which(dat.loops[,class_col_1] == cat),])
	}
	if (useNA) {
	    x = table(c(tmp[,"te_spec_l"], tmp[,"te_spec_r"]), useNA="always")
	} else {
            x = table(c(tmp[,"te_spec_l"], tmp[,"te_spec_r"]))
	}
	for (spec in c("Human", "Mouse", "Shared")) {
	    ret[i,spec] = x[spec]
	}
	if (useNA) {
	    # Add NA value counts to Shared category
	    #ret[i,"Shared"] = ret[i,"Shared"] + x[which(is.na(names(x)))]
	    ret[i,"No_Enr"] = x[which(is.na(names(x)))]
	}
    }
    ret$class = factor(ret$class, levels=cats)
    return(ret)
}

compile_class_bigwig_scores = function(dat.loops, dat.te, bigwig, cell_1, cell_2=NULL,  size=1000, cats = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")) {
    # Compile bigwig scores from the given bigwig file for TE-derived anchors in the dat.loops.
    # Return a data frame where each row is a vector of column means for windows of the given
    # size around the CTCF peak location for the TE-derived loop anchor (in dat.te)
    require(rtracklayer)    

    ret = as.data.frame(matrix(nrow=length(cats), ncol=size+1))
    colnames(ret) = c("cat", 1:size)

    class_col_1 = paste("class", cell_1, sep="_")
    if (!is.null(cell_2)) {
	class_col_2 = paste("class", cell_2, sep="_")
    }
        
    for (i in 1:length(cats)) {
    	cat = cats[i]
        x = dat.loops[which(dat.loops$te_derived & dat.loops[,class_col_1] == cat),]
	if (!is.null(cell_2)) {
	    x = rbind(x, dat.loops[which(dat.loops$te_derived & dat.loops[,class_col_2] == cat),])
	}
	for (j in 1:nrow(x)) {
	    peak = list(dat.te[which(dat.te$id_cp == x[j,"id_cp"]),])
	    for (p in peak) {
	        #pp = round(mean(c(p$chromstart_rp, p$chromend_rp)))
		pp = p$peak
	        if (!exists("ints")) {
	            ints = as.data.frame(matrix(nrow=1, ncol=3))
        	    colnames(ints) = c("chrom", "start", "end")
	    	    ints["chrom"] = dat.te[which(dat.te$id_cp == x[j,"id_cp"]), "chrom"]
	    	    ints["start"] = pp
		    ints["end"] = pp
	    	} else {
		    ints = rbind(ints, data.frame(chrom = dat.te[which(dat.te$id_cp == x[j,"id_cp"]), "chrom"], start = pp, end = pp))
	    	}
	    }
	}

	flen = round(size / 2)
	# Get data from the bigWig file
	dat = import.bw(con = BigWigFile(bigwig),
                        selection = flank(makeGRangesFromDataFrame(ints, ignore.strand=TRUE),
                                          flen, both=TRUE),
		        as = "NumericList")
	# Convert bigWig data to standard dataframe form
	dat = lapply(dat, drop)
	dat = rbind.data.frame(lapply(dat, unlist))
	dat = apply(dat, 1, unlist)
	rownames(dat) = unlist(lapply(rownames(dat), function(x){unlist(strsplit(x, "[.]"))[2]}))

	# Calculate col means and store the output vector in ret
	ret[i,2:ncol(ret)] = colMeans(dat)
	ret[i,"cat"] = cat
    }
    return(ret)
}

convert_to_long = function(dat, cats = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")) {
    # Convert wide-format matrix of bigwig values to long-format.
    for (i in 1:length(cats)) {
        cat = cats[i]
	x = dat[which(dat$cat == cat),2:ncol(dat)]
	tmp = data.frame(cat = cat, position = 1:ncol(x), value = unlist(x[1,]))
	if (!exists("ret")) {
	    ret = tmp
	} else {
	    ret = rbind(ret, tmp)
	}
    }
    return(ret)
}

wilcox_age_matrix = function(dat, cell_1, cell_2=NULL, alt="t", cats = c("C", "B2", "B1", "B0", "N1A", "N1B", "N0")) {
    # Build a matrix of Wilcoxon p-values comparing estimated TE ages between all
    # pairs of conservation classes present in a dataset.

    class_col_1 = paste("class", cell_1, sep="_")
    if (!is.null(cell_2)) {
        class_col_2 = paste("class", cell_2, sep="_")
    }

    cats = factor(cats, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
    cats = cats[order(cats)]
    ret = matrix(ncol=length(cats), nrow=length(cats))
    # Only fill in the upper triangle of the matrix
    for (i in 1:(length(cats)-1)) {
        for (j in (i+1):length(cats)) {	
            class_1 = cats[i]
            class_2 = cats[j]
	    scores_1 = dat[which(dat[,class_col_1] == class_1), "estAge"]
            scores_2 = dat[which(dat[,class_col_1] == class_2), "estAge"]
	    if (!is.null(cell_2)) {
    	        scores_1 = dat[which(dat[,class_col_1] == class_1 | dat[,class_col_2] == class_1), "estAge"]
            	scores_2 = dat[which(dat[,class_col_1] == class_2 | dat[,class_col_2] == class_2), "estAge"]
	    } else {
	        scores_1 = dat[which(dat[,class_col_1] == class_1), "estAge"]
            	scores_2 = dat[which(dat[,class_col_1] == class_2), "estAge"]
	    }
            # By default, uses alternative "two.sided" -- just testing for inequality
            f = wilcox.test(scores_1, scores_2, alternative=alt)
            ret[i,j] = f$p.value
        }
    }    
    colnames(ret) = rownames(ret) = cats
    return(ret)
}
