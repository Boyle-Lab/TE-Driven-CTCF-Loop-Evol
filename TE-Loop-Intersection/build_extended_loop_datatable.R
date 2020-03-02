# R functions for building the data tables.

add_te_intersections = function(dat.loops, dat.te) {
    # Adds TE intersection columns to loop data
    # Loop data is assumed to have a single row
    # for each annotated ChIA-pet loop, with colnames
    # c("id_cp", "chrom1", "start1", "end1", "start2", "end2")
    # Assumes that dat.te and dat.loops only contain data
    # for matching cell types!

    dat = dat.loops
    dat$te_right = dat$te_left = FALSE

    for (id_cp in dat$id_cp) {    	
    	x = dat.te[which(dat.te$id_cp == id_cp),]
	if (nrow(x) > 0) {
	    loop = dat[which(dat$id_cp == id_cp),]
	    # Compare the CTCF peak locations to the loop start/end coordinates to figure out which anchor a TE insertion represents.
	    # This just sets the boolean value for the appropriate indicator column to true!
	    for (row in 1:nrow(x)) {
	    	#print(x[row,])
		peak = x[row,"peak_ctcf"]
		if (peak >= loop$start1 &&
		    peak <= loop$end1) {
		    # Left anchor
		    dat[which(dat$id_cp == id_cp), "te_left"] = TRUE
		} else if (peak >= loop$start2 &&
                           peak <= loop$end2) {
	            # Right anchor
		    dat[which(dat$id_cp == id_cp), "te_right"] = TRUE
	        }
	    }
	}
    }
    return(dat)
}

add_motif_orientations = function(dat.loops, dat.motifs) {
    # Adds data regarding motif orienations within loop anchors.
    # NA = no motif, + = plus-strand motif, - = minus strand motif.
    # Assumes a single row in dat.loops for each ChIA-pet loop and
    # that both data frames contain data for one matching cell type.

    dat = dat.loops
    dat$motif_right = dat$motif_left = NA

    for	(id_cp in dat$id_cp) {
        x = dat.motifs[which(dat.motifs$id_cp == id_cp),]
	if (nrow(x) > 0) {
            loop = dat[which(dat$id_cp == id_cp),]
	    # Compare the midpoint of the motif to the start and end coordinates of the loop anchors to determine which end the motif resides in.
	    # If we find ANY + motif in the left anchor, we call a + motif for that anchor (even if it isn't the strongest motif).
	    # Likewise, if we fine ANY - motif in the right anchor, we call a - motif there.
	    for (row in 1:nrow(x)) {
            	#print(x[row,])
                peak = mean(x[row,"start_mt"], x[row,"end_mt"])
		if (peak >= loop$start1 &&
                    peak <= loop$end1) {
                    # Left anchor
		    if (!is.na(dat[which(dat$id_cp == id_cp), "motif_left"]) &&
		        dat[which(dat$id_cp == id_cp), "motif_left"] == "+") {
			next
		    } else {
                        dat[which(dat$id_cp == id_cp), "motif_left"] = x[row, "strand_mt"]
		    }
                } else if (peak >= loop$start2 &&
                           peak <= loop$end2) {
                    # Right anchor
		    if (!is.na(dat[which(dat$id_cp == id_cp), "motif_right"]) &&
                        dat[which(dat$id_cp == id_cp), "motif_right"] == "-") {
			next
		    } else {
                        dat[which(dat$id_cp == id_cp), "motif_right"] = x[row, "strand_mt"]
		    }
                }
	    }
	}	    
    }
    return(dat)
}

get_motif_data_table = function(dat) {
    # Prepare a table of summary data on motif layouts found in subsets of the loop data.
    ret = as.data.frame(matrix(ncol = 6, nrow = 3))
    colnames(ret) = c("Convergent", "Tandem", "Divergent", "Possible Convergent", "Possible Divergent", "No Data")
    rownames(ret) = c("Two TE", "One TE", "No TE")
    
    tmp = dat[which(dat$te_left & dat$te_right),]
    ret["Two TE",] = c(nrow(tmp[which(tmp$motif_left == "+" & tmp$motif_right == "-"),]),
    	               nrow(tmp[which(tmp$motif_left == tmp$motif_right) ,]),
		       nrow(tmp[which(tmp$motif_left == "-" & tmp$motif_right == "+"),]),
		       nrow(tmp[which( (tmp$motif_left == "+" & is.na(tmp$motif_right)) | (is.na(tmp$motif_left) & tmp$motif_right == "-") ),]),
		       nrow(tmp[which( (tmp$motif_left == "-" & is.na(tmp$motif_right)) | (is.na(tmp$motif_left) & tmp$motif_right == "+") ),]),
		       nrow(tmp[which(is.na(tmp$motif_left) & is.na(tmp$motif_right)),]))
		       
    tmp	= dat[which( (dat$te_left | dat$te_right) & !(dat$te_left & dat$te_right) ),]
    ret["One TE",] = c(nrow(tmp[which(tmp$motif_left == "+" & tmp$motif_right == "-"),]),
              	       nrow(tmp[which(tmp$motif_left == tmp$motif_right) ,]),
                       nrow(tmp[which(tmp$motif_left == "-" & tmp$motif_right == "+"),]),
                       nrow(tmp[which( (tmp$motif_left == "+" & is.na(tmp$motif_right)) | (is.na(tmp$motif_left) & tmp$motif_right == "-") ),]),
                       nrow(tmp[which( (tmp$motif_left == "-" & is.na(tmp$motif_right)) | (is.na(tmp$motif_left) & tmp$motif_right == "+") ),]),
                       nrow(tmp[which(is.na(tmp$motif_left) & is.na(tmp$motif_right)),]))
		       
    tmp = dat[which(!dat$te_left & !dat$te_right),]
    ret["No TE",] = c(nrow(tmp[which(tmp$motif_left == "+" & tmp$motif_right == "-"),]),
		       nrow(tmp[which(tmp$motif_left == tmp$motif_right) ,]),
                       nrow(tmp[which(tmp$motif_left == "-" & tmp$motif_right == "+"),]),
		       nrow(tmp[which( (tmp$motif_left == "+" & is.na(tmp$motif_right)) | (is.na(tmp$motif_left) & tmp$motif_right == "-") ),]),
		       nrow(tmp[which( (tmp$motif_left == "-" & is.na(tmp$motif_right)) | (is.na(tmp$motif_left) & tmp$motif_right == "+") ),]),
		       nrow(tmp[which(is.na(tmp$motif_left) & is.na(tmp$motif_right)),]))
		       
    return(ret)
}
