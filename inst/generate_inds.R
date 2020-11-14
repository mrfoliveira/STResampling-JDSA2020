library(doParallel)

# CHANGE NUMBER OF CORES
NCORES <- 10
NUM_SPLITS <- 10
# note that NCORES x NUM_THREADS will be actually used when running ranger
cat(paste("\nUsing", NCORES, "cores and up to", NCORES, "ranger threads\n\n"))
registerDoParallel(cores=NCORES)

#--------------------------------------

# PATHS

DATA_PATH <-  "../data/"
UTILS_PATH <- "R/"

#--------------------------------------

# LOADING PACKAGE CODE

install.packages("STResamplingExt_0.5-2.tar.gz", repos=NULL)

if(!("STResamplingExt") %in% installed.packages()){
  tosource <- list.files(UTILS_PATH, full.names = TRUE)
  for(f in tosource) source(f)
}else{
  library(STResamplingExt)
}

#--------------------------------------

cat("Loading data sets...\n")
load(paste0(DATA_PATH, "dfs_fixed.Rdata"))
# load(paste0(DATA_PATH, "20200522_inds_df.Rdata"))
dfnms <- c("MESApol", "NCDCPprec", "TCEQOozone",
           "TCEQTtemp", "TCEQWwind", "RURALpm10",
           "BEIJno", "BEIJpm10", "BEIJwind", "BEIJpm25")

inds_df <- list()
for(i in 1:length(data_list)){
  dfnm <- dfnms[i]
  
  ALPHA <- 0.25
  BETAS <- c(0.0250, 0.0375, 0.0500)
  
  if(!(dfnm %in% names(inds_df))){
    
    cat(paste("\n", dfnm))

    cat("Get spatio-temporal indicators...\n")
    ind_df <- get_full_indicators(data_list[[dfnm]]$df, data_list[[dfnm]]$stations,
                                  k=8, var="value",
                                  betas=BETAS, alpha=ALPHA,
                                  stats = c("mean", "weighted.mean", "sd"), 
                                  ratios2add = c(TRUE,TRUE,FALSE),
                                  parallel=TRUE, nsplits=NUM_SPLITS,
                                  time_id="time", site_id="station") 
    
    # fix formats
    ind_df <- as.data.frame(ind_df)
    if(grepl("BEIJ", dfnm)) ind_df$time <- lubridate::ymd_hms(ind_df$time)
    else ind_df$time <- lubridate::ymd(ind_df$time)
    
    ind_df <- list(df=ind_df, alpha=ALPHA, betas=BETAS, stations=data_list[[dfnm]]$stations)
    
    cat("\nSaving indicator data...\n")
    if(!dir.exists(paste0(DATA_PATH, "20200522/")))
      dir.create(paste0(DATA_PATH, "20200522/"))
      
    save(ind_df, file=paste0(DATA_PATH, "20200522/20200522_inds_df", dfnm,".Rdata"))
    inds_df[[dfnm]] <- ind_df 
  }
}
save(inds_df, file=paste0(DATA_PATH, "20200522_inds_df.Rdata"))

