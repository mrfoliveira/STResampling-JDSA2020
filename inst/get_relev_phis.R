d_labels <- c(10, 20, 30, 31, 32, 40, 50, 51, 52, 53)
d_names <- c("MESApol", "NCDCPprec", "TCEQOozone",
             "TCEQTtemp", "TCEQWwind", "RURALpm10",
             "BEIJno", "BEIJpm10", "BEIJwind", "BEIJpm25")

relevs <- list()
for(d in 1:10){
  print(d)
  
	load(paste0("../results/20200522/20200522_res_internal_ext_df", d,"_10repeats.Rdata"))
	
	rel <- res$results$rep1[[1]]$rawRes[[1]]$call.$internal.eval.pars$rel
	thr <- res$results$rep1[[1]]$rawRes[[1]]$call.$internal.eval.pars$thr
	cf <- res$results$rep1[[1]]$rawRes[[1]]$call.$internal.eval.pars$cf
	nORp <- res$results$rep1[[1]]$rawRes[[1]]$call.$nORp

	relevs[[d]] <- list()
	for(f in 1:length(res$results$rep1[[1]]$rawRes)){
	  
	  relevs[[d]][[f]] <- list()
	  for(type in c("rm.na", "with.na")){
	    
	    if(type == "with.na"){
	      y_train <- res$results$rep1[[1]]$rawRes[[f]]$train_y    
	    }else{
	      train_df <- inds_df[[d_names[d]]]$df[res$results$rep1[[1]]$rawRes[[f]]$train_idxs, ]
	      old_size <- nrow(train_df)
	      train_df <- centralImputTrNAs(train_df, nORp)
	      new_size <- nrow(train_df)
	      y_train <- train_df[ , "value", drop=T]
	    }
	    
	    y <- res$results$rep1[[1]]$rawRes[[f]]$results$trues
	    
	    ph <- STResamplingExt:::get_phi_control(y = y_train, rel = rel, cf = cf)
	    ls <- uba::loss.control(y_train)
	    
	    phi.trues <- UBL::phi(y,control.parms = ph)
	    phi.train <- UBL::phi(y_train,control.parms = ph)
	    
	    relevs[[d]][[f]][[type]] <- data.frame(num.relev = length(which(phi.trues>=thr)), perc.relev = length(which(phi.trues>=thr)) / length(y), 
	                                   num.relev.train = length(which(phi.train>=thr)), perc.relev.train = length(which(phi.train>=thr)) / length(y_train),
	                                   num.total.train = length(y_train))	
	    
	  }
	  relevs[[d]][[f]] <- dplyr::bind_rows(relevs[[d]][[f]], .id="pre.process")
	}
	relevs[[d]] <- dplyr::bind_rows(relevs[[d]], .id="fold")
}
relevs <- dplyr::bind_rows(relevs, .id="data") %>%
	dplyr::mutate(data = as.factor(d_labels[as.numeric(data)]),
		fold = as.numeric(fold)) %>%
  filter(type == "with.na")

save(relevs, file="../results/fold_relevs_2.Rdata")

