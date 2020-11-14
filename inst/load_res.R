
sumRes <- list()
for(d in 1:10){
  load(paste0("../results/20200227_res_internal_ext_df", d, ".Rdata"))
  
  sumRes[[d]] <- lapply(res$results, 
                   function(x){ 
                     ires <- as.data.frame(x$evalRes);
                     colnames(ires)[which(colnames(ires)=="bias")] <- "rmse.bias";
                     colnames(ires)[which(colnames(ires)=="var")] <- "rmse.var";
                     
                     chosen <- dplyr::bind_rows(lapply(x$rawRes, 
                                  function(y) as.data.frame(y$resample_chosen)));
                     colnames(chosen) <- paste0("chosen.", colnames(chosen));
                     idf <- cbind(cbind(data.frame(df=d, fold=1:9), ires), chosen);
                     idf
                     })
  sumRes[[d]] <- dplyr::bind_rows(lapply(1:length(sumRes[[d]]), 
                                    function(i) 
                                      cbind(sumRes[[d]][[i]], res$fixed.grid[i,,drop=F])))
}
sumRes <- dplyr::bind_rows(sumRes)

save(sumRes, file="../results/20200227_sumRes_internal_ext.Rdata")



sumRes <- list()
for(d in 1:10){
  load(paste0("../results/20200227_res_oracle_ext_df", d, ".Rdata"))
  
  sumRes[[d]] <- lapply(res$results, 
                        function(x){ 
                          ires <- as.data.frame(x$evalRes);
                          colnames(ires)[which(colnames(ires)=="bias")] <- "rmse.bias";
                          colnames(ires)[which(colnames(ires)=="var")] <- "rmse.var";
                          
                          idf <- cbind(data.frame(df=d, fold=1:9), ires);
                          idf
                        })
  
  sumRes[[d]] <- dplyr::bind_rows(lapply(1:length(sumRes[[d]]), 
                                         function(i) 
                                           cbind(sumRes[[d]][[i]], res$oracle.grid[i,,drop=F])))
}
sumRes <- dplyr::bind_rows(sumRes)

save(sumRes, file="../results/20200227_sumRes_oracle_ext_marineye.Rdata")


#for(d in 1:10){
# cat(d)
# # checking that data sets match
# print(all(inds_df[[d]]$df[res$results$rep1[[1]]$rawRes[[1]]$test_idxs,"value"] == res$results$rep1[[1]]$rawRes[[1]]$results$trues))
#
# for(j in 1:length(res$results)){
#	for(k in 1:length(res$results[[j]])){
#       # saving memory space   
#		for(ii in 1:length(res$results[[j]][[k]]$rawRes)){
# 			res$results[[j]][[k]]$rawRes[[ii]]$internal_call.$train <- NULL	
#			res$results[[j]][[k]]$rawRes[[ii]]$internal_call.$test <- NULL	
#		}
# 	}
# }
# cat("and saving ... \n")
# save(res, file=paste0("../results/20200408_res_internal_ext_df", d, "_10repeats.Rdata"))
#}
  	

# loading internal validation results with repetition

sumRes <- list()
for(d in 1:10){
  cat(d)
  load(paste0("../results/20200522/20200522_res_internal_ext_df", d, "_10repeats.Rdata"))
  cat(" loaded ... ")

  sumRes[[d]] <- list()
  for(rep in 1:length(res$results)){
  	#sumRes[[d]][[rep]] <- lapply(res$results[[rep]], 
    #                    function(x){ 
    sumRes[[d]][[rep]] <- list()
    for(exp in 1:length(res$results[[rep]])){                          
    					  x <- res$results[[rep]][[exp]]
                          if( is.list(x) ){
                        		ires <- as.data.frame(x$evalRes);
		                        colnames(ires)[which(colnames(ires)=="bias")] <- "rmse.bias";
		                        colnames(ires)[which(colnames(ires)=="var")] <- "rmse.var";
		                          
		                        chosen <- dplyr::bind_rows(lapply(x$rawRes, 
		                                                            function(y) as.data.frame(y$resample_chosen)));
		                        colnames(chosen) <- paste0("chosen.", colnames(chosen));
		                        idf <- cbind(cbind(cbind(data.frame(df=d, fold=1:9, rep=rep), ires), 
                              chosen), res$fixed.grid[exp,,drop=F]);
	                      }else{
	                      	idf <- NULL
	                      }
                          sumRes[[d]][[rep]][[exp]] <- idf
                          #idf
                        #})
                        }
    sumRes[[d]][[rep]] <- dplyr::bind_rows(sumRes[[d]][[rep]])
  }
  sumRes[[d]] <- dplyr::bind_rows(sumRes[[d]])
}
sumRes <- dplyr::bind_rows(sumRes)

save(sumRes, file="../results/20200522_sumRes_internal_ext_10repeats_nitro.Rdata")


# loading oracle results with repetition

sumRes <- list()
for(d in 1:10){
  cat(paste("Data", d))
  load(paste0("../results/20200522/20200522_res_oracle_ext_df", d, "_10repeats_ratio.Rdata"))
  cat(" loaded... ")

  sumRes[[d]] <- list()
  for(rep in 1:length(res$results)){
  	cat(paste0(rep, " "))
  	sumRes[[d]][[rep]] <- list()
  	for(exp in 1:length(res$results[[rep]])){
    #sumRes[[d]][[rep]] <- lapply(res$results[[rep]], 
    #                    function(x){
    						# cat(exp)
    						x <- res$results[[rep]][[exp]]
	                        if(is.list(x)){ 
	                          ires <- as.data.frame(x$evalRes);
	                          colnames(ires)[which(colnames(ires)=="bias")] <- "rmse.bias";
	                          colnames(ires)[which(colnames(ires)=="var")] <- "rmse.var";
	                          
	                          idf <- cbind(dplyr::bind_cols(list(data.frame(df=rep(d,9), fold=1:9, rep=rep(rep, 9)), 
	                          	ires)), as.data.frame(res$oracle.grid[exp,,drop=F]));
	                      	}else{
	                      	  idf <- NULL	
	                      	}
	                      	sumRes[[d]][[rep]][[exp]] <- idf
                          #idf
                        #})

  	}
    sumRes[[d]][[rep]] <- dplyr::bind_rows(sumRes[[d]][[rep]])                                         
  }
  sumRes[[d]] <- dplyr::bind_rows(sumRes[[d]])
}
sumRes <- dplyr::bind_rows(sumRes)

save(sumRes, file="../results/20200522_sumRes_oracle_ext_10repeats_ratio_marineye.Rdata")
