library(ranger)
library(earth)
library(rpart)

#--------------------------------------

# PATHS

DATA_PATH <- system.file("inst/extdata/", package = "STResamplingJDSA")
UTILS_PATH <- "R" # if package not installed
RESULTS_PATH <- "inst/"
if(!dir.exists(RESULTS_PATH)) dir.create(RESULTS_PATH)

#--------------------------------------

# LOADING PACKAGE CODE

if(!("STResamplingExt") %in% installed.packages()){
  devtools::install_github("mrfoliveira/STResampling-JDSA2020",ref="master")
  library(STResamplingJDSA)
  
  # tosource <- list.files(UTILS_PATH, full.names = TRUE)
  # for(f in tosource) source(f)
}else{
  library(STResamplingJDSA)
}

#--------------------------------------

# FIXED PARAMETERS

d_names <- c("MESApol", "NCDCPprec", "TCEQOozone",
  "TCEQTtemp", "TCEQWwind", "RURALpm10",
  "BEIJno", "BEIJpm10", "BEIJwind", "BEIJpm25")

SEED <- 1234
NREPS <- 10

THR_REL <- 0.9
REL <- "auto"
CF <- 1.5
BETA <- 1
MIN_TRAIN <- 2
NORP <- 0.2

# OTHER PARAMETERS
EST_PARS <- list(nfolds = 10, 
                 window = "growing", 
                 fold.alloc.proc = "Tblock_SPall", 
                 removeSP = FALSE, 
                 time="time", 
                 site_id="station",
                 .keepTrain = TRUE,
                 .parallel = FALSE, 
                 .verbose = FALSE)


WF_PARS_BASE <- list(min_train=MIN_TRAIN, handleNAs="centralImputNAs", nORp = NORP)
RS_PARS_BASE <- list(thr.rel=THR_REL, rel = REL, cf = CF)
EVAL_PARS <- list(eval.function = "eval_stats", rel = REL, cf = CF, 
                  thr = THR_REL, beta = BETA, .keptTrain = TRUE)

# PARAMETER VALUES TO TEST

models <- c("rpart", "earth", "ranger")

cpercs <- list()
cpercs[["under"]] <- c(0.2, 0.4, 0.6, 0.8, 0.95)
cpercs[["over"]] <- c(0.5, 1, 2, 3, 4)
cpercs[["gauss"]] <- cpercs[["over"]]
alphas <- c(0, 0.25, 0.5, 0.75, 1)
ptypes <- c("orig", "addPhiLook", "addNoPhi")
perts <- 0.1
epsilons <- 1E-4


#--------------------------------------

# BUILD PARAMETER GRID

grids <- list()
for(rsfun in c("under", "over", "gauss")){
  repl <- ifelse(rsfun == "under", FALSE, TRUE)
  
  base.params <- list(type = rsfun, model = models, repl = repl)
  if(rsfun == "gauss") base.params <- c(base.params, list(pert=perts))
  bias.params <- list(ptype = ptypes, alpha = alphas, epsilon = epsilons, bias = TRUE)
  
  grids <- dplyr::bind_rows(grids, 
                            do.call("expand.grid", c(base.params, 
                                                     list(C.perc = cpercs[[rsfun]], 
                                                          arg.type = "C.perc", bias=FALSE),
                                                     stringsAsFactors = FALSE)),
                            do.call("expand.grid", c(base.params, bias.params, 
                                                     list(C.perc = cpercs[[rsfun]],
                                                          arg.type = "C.perc"), 
                                                     stringsAsFactors = FALSE)))
}
# join grids
param.grid <- dplyr::bind_rows(grids)
# order rows
param.grid <- param.grid[order(param.grid$model, param.grid$type, param.grid$bias), ]
param.grid <- as.data.frame(param.grid, stringsAsFactors=FALSE)
