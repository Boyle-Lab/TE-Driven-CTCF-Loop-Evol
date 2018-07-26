build_wilcox_matrix = function(dat, cell_t, alt="t") {
    # Build a matrix of Wilcoxon p-values comparing scores between all
    # pairs of conservation classes present in a dataset. The diagonal
    # of the matrix contains comparisons between TE-derived and non-TE-derived
    # fractions of the given class.
    class_col = paste("class", cell_t, sep='_')
    classes = unique(dat[,class_col])
    classes = factor(classes, levels=c("C", "B2", "B1", "B0", "N1A", "N1B", "N0"))
    classes = classes[order(classes)]
    ret = matrix(ncol=length(classes), nrow=length(classes))
    # Only fill in the upper triangle of the matrix
    for (i in 1:length(classes)) {
        for (j in i:length(classes)) {
	    class_1 = classes[i]
	    class_2 = classes[j]
	    scores_1 = dat[which(dat[,class_col] == class_1), "IAB"]
	    scores_2 = dat[which(dat[,class_col] == class_2), "IAB"]
	    if (i == j) {
	        # Diagonal of the matrix. Test TE-derived versus non-TE-derived.
		scores_1 = dat[which(dat[,class_col] == class_1 & dat$te_derived == FALSE), "IAB"]
		scores_2 = dat[which(dat[,class_col] == class_1 & dat$te_derived == TRUE), "IAB"]
	    }
	    # By default, uses alternative "two.sided" -- just testing for inequality
	    f = wilcox.test(scores_1, scores_2, alternative=alt)
	    ret[i,j] = f$p.value
	}
    }    
    colnames(ret) = rownames(ret) = classes
    return(ret)
}

