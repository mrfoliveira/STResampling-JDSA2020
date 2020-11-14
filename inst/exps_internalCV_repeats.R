
source(system.file("inst/set_parameters_general.R", package = "STResamplingJDSA"))
source(system.file("ints/set_parameters_internal.R", package = "STResamplingJDSA"))


tuning.grid <- tuning.grid[c(which(tuning.grid$model!="ranger"), 
	which(tuning.grid$model=="ranger")), ]

# RUNNING EXPERIMENTS

NCORES <- 24
for(d in 1:10){
    
  # CREATE CLUSTER
  if(.PARALLEL){
    cl <- makeCluster(NCORES)
    registerDoParallel(cores=NCORES)   
  }
  
  # SET SEED
  set.seed(SEED)
  
  # LOAD DATA
  dfnm <- d_names[d]
  cat(paste("\n----- DATA SET:", dfnm, "\n\n"))
  load(paste0(DATA_PATH, "inds_df", dfnm,".Rdata"))
  stations <- ind_df$stations
  ind_df <- as.data.frame(ind_df$df)
  
  # MULTIPLE REPITITIONS
  res <- foreach::foreach(s = 1:NREPS, .errorhandling="stop") %do% {
    
    cat(paste("\nREPEAT:", s, "...\n\n"))
    
    # CYCLE THROUGH DIFFERENT PARAMETRIZATIONS
    results <- foreach::foreach(i = 1:nrow(tuning.grid), .errorhandling="pass") %dopar% {
		
		cat(paste(i, " "))
      # these are the fixed parameters
      common.params <- as.list(tuning.grid[i,])[which(!is.na(tuning.grid[i,]))]
      this.arg.type <- common.params$arg.type
      common.params$arg.type <- NULL
      # assuming the rest of the parameters are RS_PARS
      rmCols <- which(names(common.params) %in% c("model"))
      
      # cat(paste0("\nTuning ", i, ": ", 
      #           paste(names(common.params), common.params, sep=" = ", collapse = "; ")))
      
      # add parameters specific to random forest
      if(common.params$model == "ranger"){
        WF_PARS <-   c(INT_WF_PARS_BASE, list(num.trees = 250, 
                                              num.threads=NUM_THREADS, verbose=FALSE))
      }else{
        WF_PARS <- INT_WF_PARS_BASE
      } 
      # add parameters specific to biased resampling
      RS_PARS <- c(RS_PARS_BASE, common.params[-rmCols])
      if(!is.null(common.params$bias)){
        if(common.params$bias) RS_PARS <-  c(RS_PARS, list(sites_sf = stations)) 
      } 
      
      # now we find what are the parameters to tune
      this.params <- param.grid[which(param.grid$arg.type == tuning.grid[i, "arg.type"]),
                                -which(colnames(param.grid) == "arg.type")]
      # get lines that correspond to the fixed parameters
      for(j in 1:length(common.params))
        this.params <- this.params[which(this.params[[names(common.params)[j]]] == common.params[j]), ]
      # remove columns that are all just NA
      resample.grid <- this.params[ , which(colSums(is.na(this.params))<nrow(this.params))]
      # remove column names that define what groups have separate tuning
      resample.grid <- resample.grid[ , -which(colnames(resample.grid) %in% colnames(tuning.grid))]
      
      if(common.params$model == "ranger" & s > 1){
        # not repeating RF experiments more than once
        res <- "skip"
      }else{
        # get results with these parameters
        res <- estimates(data = ind_df, form = value~., 
                         estimator = "prequential_eval", est.pars = EST_PARS, 
                         workflow = "internal_workflow", 
                         wf.pars = c(list(model = common.params[["model"]],
                                          resample.grid = resample.grid,
                                          resample.pars = RS_PARS),
                                     WF_PARS), 
                         evaluator = "evaluate", 
                         eval.pars = EVAL_PARS, 
                         seed = NULL)
      }
      
     res

    }

    results

  }

  names(res) <- paste0("rep", 1:length(res))
  res <- list(results = res, fixed.grid = tuning.grid) 
  
  save(res, file = paste0(RESULTS_PATH,"res_internal_ext_df", d, "_10repeats.Rdata"))

  
  if(.PARALLEL) stopCluster(cl)
 
}
