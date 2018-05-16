calc_rpt_enrichments = function(dat.rpt, dat.words, dat.bg, species) {
    ret = as.data.frame(matrix(nrow=length(unique(dat.rpt$rpt)), ncol=6))
    colnames(ret) = c("count_fg", "specific_fg", "count_bg", "specific_bg", "pval", "pval_cor")
    rownames(ret) = unique(dat.rpt$rpt)

    for (rpt in rownames(ret)) {
        ret[rpt,"count_fg"] = nrow(dat.rpt[which(dat.rpt$rpt == rpt),])
	ret[rpt,"count_bg"] = nrow(dat.bg) - ret[rpt,"count_fg"]
    	if (species == "hg19") {
            ret[rpt,"specific_fg"] = nrow(dat.rpt[which(dat.rpt$rpt == rpt & dat.rpt$word %in% rownames(dat.words[which(dat.words$nocc.hg19 >= 8 & dat.words$orr >= 2),])),])
	    ret[rpt,"specific_bg"] = nrow(dat.bg[which(dat.bg$word %in% rownames(dat.words[which(dat.words$nocc.hg19 >= 8 & dat.words$orr >= 2),])),]) - ret[rpt,"specific_fg"]
	} else {
            ret[rpt,"specific_fg"] = nrow(dat.rpt[which(dat.rpt$rpt == rpt & dat.rpt$word %in% rownames(dat.words[which(dat.words$nocc.mm9 >= 8 & dat.words$orr <= -2),])),])
	    ret[rpt,"specific_bg"] = nrow(dat.bg[which(dat.bg$word %in% rownames(dat.words[which(dat.words$nocc.mm9 >= 8 & dat.words$orr <= -2),])),]) - ret[rpt,"specific_fg"]
	}
        FG1 = ret[rpt,"specific_fg"]
        FG2 = ret[rpt,"count_fg"] - ret[rpt,"specific_fg"]
        BG1 = ret[rpt,"specific_bg"]
        BG2 = ret[rpt,"count_bg"] - ret[rpt,"specific_bg"]
        tab =  t(matrix(c(FG1, FG2, BG1, BG2), nrow=2, ncol=2))
	pv = fisher.test(tab, alternative = "g")
        ret[rpt, "pval"] = pv$p.value
    }
    
    ret$pval_cor = p.adjust(ret$pval, method="bonferroni")
    return(ret)
}
