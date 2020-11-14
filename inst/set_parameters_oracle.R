# test-run
# param.grid <- param.grid[which(param.grid$model=="rpart" & param.grid$alpha %in% c(NA, 0.5) 
#                                & param.grid$C.perc %in% c(2,0.6,NA),]
param.grid <- dplyr::bind_rows(param.grid, 
                               as.data.frame(expand.grid(type="none", model=models),stringsAsFactors=FALSE))

# forget about RF
param.grid <- param.grid[ which(param.grid$arg.type %in% c("C.perc", NA)), 
                          -which(colnames(param.grid) %in% c("arg.type")) ]


.PARALLEL <- TRUE

if(.PARALLEL){
  # PARALLELIZATION
  NCORES <- 24
  NUM_SPLITS <- NCORES
  NUM_THREADS <- 1
  library(doParallel)
  cat(paste("\nUsing", NCORES, "cores and up to", NCORES*NUM_THREADS, "ranger threads\n\n"))
  # registerDoParallel(cores=NCORES)
}