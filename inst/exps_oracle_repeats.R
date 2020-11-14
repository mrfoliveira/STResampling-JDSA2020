
source(system.file("inst/set_parameters_general.R", package = "STResamplingJDSA"))
source(system.file("ints/set_parameters_oracle.R", package = "STResamplingJDSA"))

# RUNNING EXPERIMENTS

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
  res <- foreach::foreach(s = 1:NREPS, .errorhandling="pass") %do% {
    
    cat(paste("\nREPEAT:", s, "...\n\n"))
    
    # CYCLE THROUGH DIFFERENT PARAMETRIZATIONS
    results <- foreach::foreach(i = 1:nrow(param.grid), .errorhandling="pass") %dopar% {
      
      params <- as.list(param.grid[i,])[which(!is.na(param.grid[i,]))]
      rmCols <- which(names(params) %in% c("model"))
      
      # add parameters specific to random forest
      if(params$model == "ranger"){
        WF_PARS <-   c(WF_PARS_BASE, list(num.trees = 250, 
                                          num.threads=NUM_THREADS, verbose=FALSE))
      }else{
        WF_PARS <- WF_PARS_BASE
      } 
      # add parameters specific to gaussian and no re-sampling
      RS_PARS <- if(params$type == "none") NULL else c(RS_PARS_BASE, as.list(params[-rmCols]))
      # add parameters specific to biased resampling
      if(!is.null(params$bias)){
        if(params$bias) RS_PARS <-  c(RS_PARS, list(sites_sf = stations)) 
      } 
      
      if(params$model == "ranger" & s > 1){
        # not repeating RF experiments more than once
        res <- "skip"
      }else{
        # get results with these parameters
        res <- estimates(data = ind_df, form = value~., 
                         estimator = "prequential_eval",
                         est.pars = EST_PARS, 
                         workflow = "simple_workflow", 
                         wf.pars = c(list(model = params[["model"]],
                                          resample.pars = RS_PARS),
                                     WF_PARS), 
                         evaluator = "evaluate", 
                         eval.pars = EVAL_PARS, 
                         seed = NULL)  
      }
      
      res
    }
      
    save(results, file = paste0(RESULTS_PATH,"res_oracle_ext_df", d, "_rep", s,".Rdata"))
    rm(results)
    gc()
    
    s
  }
  
  res <- list()
  for(s in 1:10){
      load(paste0(RESULTS_PATH,"res_oracle_ext_df", d, "_rep", s,".Rdata"))
    res[[s]] <- results
    file.remove(paste0(RESULTS_PATH,"res_oracle_ext_df", d, "_rep", s,".Rdata"))
  }
  res <- c(list(results = res, oracle.grid = param.grid))
  save(res, file = paste0(RESULTS_PATH,"res_oracle_ext_df", d, "_10repeats.Rdata"))
  rm(res)
  rm(results)
  gc()
  
  if(.PARALLEL) stopCluster(cl)
}
