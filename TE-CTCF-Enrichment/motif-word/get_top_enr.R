get_te_word_data = function(dat.motifs, dat.enr, te_list, species) {

    ret = as.data.frame(matrix(nrow=length(te_list), ncol=ncol(dat.enr)+3))
    colnames(ret) = c(colnames(dat.enr), "word", "total_bound_insertions", "family_wise_bound_insertions")
    rownames(ret) = te_list

    for (fam in te_list) {
        dat = summary(factor(dat.motifs[which(dat.motifs$rpt == fam),"word"]), maxsum=nrow(dat.motifs[which(dat.motifs$rpt == fam),]))
	if (length(dat) == 0) {
            next
        }
        top_word = names(dat[which(dat == max(dat, na.rm=TRUE))])
	dat2 = dat.enr[top_word,]
	if (species == "hg19") {
	    dat2 = dat2[order(dat2$orr, decreasing=TRUE),]
	} else {
	    dat2 = dat2[order(dat2$orr, decreasing=FALSE),]
	}
	n_bound = nrow(dat.motifs[which(dat.motifs$rpt == fam & dat.motifs$word == rownames(dat2[1,])),])
        ret[fam,] = c(dat2[1,], rownames(dat2[1,]), nrow(dat.motifs[which(dat.motifs$rpt == fam),]), n_bound)
    }
    return(ret)
}