do_upset_plot = function(dat.all, query_cell, target_cell, setname) {
    require(UpSetR)
    dat = dat.all[which(dat.all$cell == query_cell),]

    tcol = colnames(dat)[which(dat[1,] == target_cell)]
    t_ext = unlist(strsplit(tcol, "_"))[2]
    target_col = paste("class", t_ext, sep="_")
    maps_col_l = paste("l_maps", t_ext, sep="_")
    maps_col_r = paste("r_maps", t_ext, sep="_")

    dat_upset = c("Conserved" = nrow(dat[(dat[,target_col] == "C"),]),
              "Conserved&TE" = nrow(dat[which(dat[,target_col] == "C" & (dat$te_left | dat$te_right)),]),
              "B2" = nrow(dat[which(dat[,target_col] == "B2"),]),
              "B2&TE" = nrow(dat[which(dat[,target_col] == "B2" & (dat$te_left | dat$te_right)),]),
              "B1" = nrow(dat[which(dat[,target_col] == "B1"),]),
              "B1&TE" = nrow(dat[which(dat[,target_col] == "B1" & (dat$te_left | dat$te_right)),]),
              "B0" = nrow(dat[which(dat[,target_col] == "B0"),]),
              "B0&TE" = nrow(dat[which(dat[,target_col] == "B0" & (dat$te_left | dat$te_right)),]),
              "N1A" = nrow(dat[which(dat[,target_col] == "N1A"),]),
              "N1A&TE" = nrow(dat[which(dat[,target_col] == "N1A" & ( (dat$te_left & dat[,maps_col_l] == FALSE) | (dat$te_right & dat[,maps_col_r] == FALSE) )),]),
              "N1B" = nrow(dat[which(dat[,target_col] == "N1B"),]),
              "N1B&TE" = nrow(dat[which(dat[,target_col] == "N1B" & ( (dat$te_left & dat[,maps_col_l] == FALSE) | (dat$te_right & dat[,maps_col_r] == FALSE) )),]),
              "N0" = nrow(dat[which(dat[,target_col] == "N0"),]),
              "N0&TE" = nrow(dat[which(dat[,target_col] == "N0" & ( (dat$te_left & dat[,maps_col_l] == FALSE) | (dat$te_right & dat[,maps_col_r] == FALSE) )),]))
    return(dat_upset)
    pdf(paste("upset_", query_cell, "-to-", target_cell, ".", setname, ".pdf", sep=""))
    upset(fromExpression(dat_upset), nsets=length(dat_upset), order.by="degree", group.by="sets")
    dev.off()
}

do_bar_plot = function(dat.all, query_cell, target_cell, setname) {
    dat = dat.all[which(dat.all$cell == query_cell),]

    tcol = colnames(dat)[which(dat[1,] == target_cell)]
    t_ext = unlist(strsplit(tcol, "_"))[2]
    target_col = paste("class", t_ext, sep="_")
    maps_col_l = paste("l_maps", t_ext, sep="_")
    maps_col_r = paste("r_maps", t_ext, sep="_")

    require(ggplot2)
    require(reshape)
    dat.sets = list("Conserved" = c( nrow(dat[(dat[,target_col] == "C"),]),
               			     nrow(dat[which(dat[,target_col] == "C" & (dat$te_left | dat$te_right)),]) ),
                    "B2" = c( nrow(dat[which(dat[,target_col] == "B2"),]),
              	    	      nrow(dat[which(dat[,target_col] == "B2" & (dat$te_left | dat$te_right)),]) ),
                    "B1" = c( nrow(dat[which(dat[,target_col] == "B1"),]),
              	    	      nrow(dat[which(dat[,target_col] == "B1" & (dat$te_left | dat$te_right)),]) ),
              	    "B0" = c( nrow(dat[which(dat[,target_col] == "B0"),]),
        	    	      nrow(dat[which(dat[,target_col] == "B0" & (dat$te_left | dat$te_right)),]) ),
              	    "N1A" = c( nrow(dat[which(dat[,target_col] == "N1A"),]),
              	    	       nrow(dat[which(dat[,target_col] == "N1A" & ( (dat$te_left & dat[,maps_col_l] == FALSE) | (dat$te_right & dat[,maps_col_r] == FALSE) )),]) ),
              	    "N1B" = c( nrow(dat[which(dat[,target_col] == "N1B"),]),
		    	       nrow(dat[which(dat[,target_col] == "N1B" & ( (dat$te_left & dat[,maps_col_l] == FALSE) | (dat$te_right & dat[,maps_col_r] == FALSE) )),]) ),
              	    "N0" = c( nrow(dat[which(dat[,target_col] == "N0"),]),
		    	      nrow(dat[which(dat[,target_col] == "N0" & ( (dat$te_left & dat[,maps_col_l] == FALSE) | (dat$te_right & dat[,maps_col_r] == FALSE) )),])) )

    dat.sets = as.data.frame(t(as.data.frame(dat.sets)))
    dat.sets[,1] = dat.sets[,1] - dat.sets[,2]
    colnames(dat.sets) = c("No TE", "TE")
    dat.sets$cat = factor(rownames(dat.sets), levels = c("Conserved", "B2", "B1", "B0", "N1A", "N1B", "N0"))

    return(dat.sets)

    pdf(paste("stacked-bar_", query_cell, "-to-", target_cell, ".", setname, ".pdf", sep=""))
    ggplot(melt(dat.sets), aes(x=cat, y=value, fill=variable))+
    geom_bar(stat='identity', position='fill')
    dev.off()
}