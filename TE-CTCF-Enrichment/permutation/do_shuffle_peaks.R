require(fgpt)
require(sqldf)

shuffle_rep_data = function(cdat, rdat, nperm=10000, by_bin=TRUE) {
    # cdat = ctcf_repeats
    # rdat = genome-wide repeats
    
    fams = sort(unique(cdat$name_rmsk))
    m = matrix(nrow=length(fams), ncol=nperm)
    rownames(m) = fams
    for (i in 1:nperm) {
        cat(".")
        # Shuffle repeat names genome-wide and store the number of observations
	# for each family observed in the ctcf-repeat data
	shuf = rdat[,c("id_rmsk", "name", "tss_dist_bin")]
	#print(head(shuf))
	if (by_bin) {
            for (j in 1:6) {
	        if (length(rdat[which(rdat$tss_dist_bin == j),"name"]) == 0) {
	            next
	        }
	    	#print(j)
            	shuf[which(shuf$tss_dist_bin == j),"name"] = fyshuffle(rdat[which(rdat$tss_dist_bin == j),"name"])
            }
	} else {
	    shuf[,"name"] = fyshuffle(rdat[,"name"])
	}
	# Intesect the shuffled data with the real data
	mdat = sqldf("SELECT x.id_rmsk, x.name_rmsk AS 'name', y.name AS 'name_perm' FROM cdat x INNER JOIN shuf y ON x.id_rmsk = y.id_rmsk")
	#print(head(mdat))
	for (fam in fams) {
	   # Count the number of instances of the observed family in permuted data
    	   m[fam,i] = nrow(mdat[which(mdat$name == fam & mdat$name_perm == fam),])
	   #print(c(fam, m[fam,]))
	}
	#print(head(m))
    }
    cat("", sep="\n")
    # Store the empirical CDFs
    ret = data.frame(observed_count = rep(0, length(fams)), avg_perm_count = rep(0, length(fams)), min_perm_count = rep(0, length(fams)), max_perm_count = rep(0, length(fams)), pval = rep(1, length(fams)), pval_cor = rep(1, length(fams)), pval_bn = rep(1, length(fams)), pval_bn_cor = rep(1, length(fams)), species = rep(cdat$species[1], length(fams)), cell = rep(cdat$cell[1], length(fams)), row.names = fams, stringsAsFactors = FALSE)
    for (fam in fams) {
        n = nrow(cdat[which(cdat$name_rmsk == fam),])
	if (n > 0 && max(m[fam,]) > 0) {
	    ret[fam,"observed_count"] = n
	    ret[fam,"avg_perm_count"] = mean(m[fam,])
	    ret[fam,"min_perm_count"] =	min(m[fam,])
	    ret[fam,"max_perm_count"] =	max(m[fam,])
            f = ecdf(m[fam,])
	    ret[fam,"pval"] = 1-f(n)
	    bn = binom.test(n, nrow(cdat), mean(m[fam,]) / nrow(cdat))
	    ret[fam,"pval_bn"] = bn$p.value
        }
    }
    #print(head(ret))
    ret$pval_cor = p.adjust(ret$pval, method="bonferroni")
    ret$pval_bn_cor = p.adjust(ret$pval_bn, method="bonferroni")
    ret$type = rownames(ret)
    rownames(ret) = NULL
    return(ret)
}

do_permutations = function(cdat, rdat, nperm=10000, by_bin=TRUE) {
    i = 1
    for (species in c("hg19", "mm9")) {

        cells = c("GM12878", "K562")
        if (species == "mm9") {
	    cells = c("CH12", "MEL")
        }

        for (cell in cells) {
            print(c(species, cell))
            res = shuffle_rep_data(cdat[which(cdat$cell == cell),], rdat[which(rdat$species == species),], nperm=nperm, by_bin=by_bin)
            if (i > 1) {
                ret = rbind(ret, res)
            } else {
                ret = res
            }
	    i = i+1
        }
    }
    return(ret)
}
