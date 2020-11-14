
# ADDING INTERNAL-SPECIFIC PARAMETERS

INT_EST_PARS <- list(nfolds = 10, 
                     fold.alloc.proc = "Tblock_SPall", 
                     time="time", 
                     site_id="station",
                     .keepTrain = TRUE,
                     .parallel = FALSE,
                     .verbose = FALSE)

INT_EVAL_PARS <- EVAL_PARS
INT_EVAL_PARS$.keptTrain <- NULL

INT_WF_PARS_BASE <- c(WF_PARS_BASE,
                      list(internal.est = "kf_xval", 
                           internal.est.pars = INT_EST_PARS,
                           internal.evaluator = "int_util_evaluate", 
                           internal.eval.pars = INT_EVAL_PARS,
                           metrics = c("F1.u", "rmse_phi"), 
                           metrics.max = c(TRUE, FALSE),
                           stat="MEAN", .int_parallel = FALSE,
                           .intRes = TRUE))

GROUP_VARS <- c("type", "model", "ptype", "bias", "arg.type")

# RUNNING EXPERIMENTS

# these are the different groups that will have separate tuning
tuning.grid <- unique(param.grid[ , which(colnames(param.grid) %in% GROUP_VARS)])
rownames(tuning.grid) <- NULL

tuning.grid <- tuning.grid[ which(tuning.grid$arg.type == "C.perc"), ]


.PARALLEL <- TRUE

if(.PARALLEL){
  # PARALLELIZATION
  NCORES <- 30
  NUM_SPLITS <- NCORES
  NUM_THREADS <- 2
  library(doParallel)
  cat(paste("\nUsing", NCORES, "cores and up to", NCORES*NUM_THREADS, "ranger threads\n\n"))
  # workers <- makeCluster(NCORES)
  # registerDoParallel(cores=NCORES) #workers)
}
