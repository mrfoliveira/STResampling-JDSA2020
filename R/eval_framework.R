############################
#### EVALUATION METRICS CALCULATOR
############################


#' Calculate regression evaluation metrics
#' 
#' Calculate regression metrics for imbalanced domains.
#' @author Nuno Moniz
#' @param trues a vector of true values
#' @param preds a vector of predicted values
#' @param y_train a vector of training values
#' @param rel The relevance function which can be automatically ("auto") 
#' determined (the default) or may be provided by the user through a matrix 
#' with interpolating points.
#' @param cf phi.control coef (needed if rel = 'auto')
#' @param thr relevance threshold
#' @param beta beta in F-measure
#' @return a named vector of metrics RMSE, relevance-aware 
#' RMSE, and utility-based precision, recall and F-measure
#' 
#' @export
eval_stats <- function(trues, preds, y_train = NULL, thr=0.9, 
                       beta=1, rel = "auto", cf = 1.5) {
  
  if(all(is.na(preds)) | length(preds)==0){
    rmse = rmse_phi = prec.u = rec.u = F1.u = bias.sq = variance = NA
  }else{
   
    
    ph <- get_phi_control(y = y_train, rel = rel, cf = cf)
    ls <- uba::loss.control(y_train)
    
    u_new <- uba::util(preds,trues, ph, ls,
                       uba::util.control(umetric="P",event.thr=thr), return.uv = TRUE)
    
    phi.trues <- UBL::phi(trues,control.parms = ph)
    phi.preds <- UBL::phi(preds,control.parms = ph)
    
    pr_frm <- data.frame(Utility=u_new)
    pr_frm["phiTrues"] <- phi.trues
    pr_frm["phiPreds"] <- phi.preds
    
    rmse <- sqrt(mean((trues-preds)^2))
    rmse_phi= sqrt(mean(phi.trues[phi.trues>thr]*(trues[phi.trues>thr]-preds[phi.trues>thr])^2))
    
    bias.sq <- mean(mean(preds) - trues)^2
    variance <- mean(mean(preds) - preds)^2
    
    prec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>thr,]$phiPreds)
    rec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>thr,]$phiTrues)
    F1.u <- (1+beta) * prec.u * rec.u / ( beta^2 * prec.u + rec.u)
   
    # clean results
    if(is.na(rmse_phi)) prec.u=rec.u=F1.u=NA
    if(!is.na(rmse_phi)){
      if(is.na(prec.u)) prec.u <- 0
      if(is.na(rec.u)) rec.u <- 0
      if(is.na(F1.u)) F1.u <- 0
    }
  }
    
  c(
    rmse=rmse, rmse_phi=rmse_phi, prec.u=prec.u, rec.u=rec.u, F1.u=F1.u, rmse.bias=bias.sq, rmse.var=variance
  )
}

############################
#### EVALUATION FRAMEWORK
############################


#' Evalute the results of a predictive workflow
#' 
#' Calculate evaluation metrics from the raw results of a workflow
#' @param wfRes a data frame (or list of data frames) containing the results of
#' a predictive workflow with columns \code{trues} and \code{preds} containing
#' the real and predicted values, respectively
#' @param eval.function the function to be used to calculate error metrics from \code{wfRes}
#' @param .keptTrain a Boolean indicating whether \code{.keepTrain} was
#' set to \code{TRUE} in calls to estimation methods. Only useful
#' if evaluation metrics need training data.
#' @param ... parameters to pass to \code{eval.function}
#'
#' @return The results (or a list of results) of \code{eval.function} applied to 
#' the data frame (or list of data frames) in \code{wfRes}
#' 
#' @export
evaluate <- function(wfRes, eval.function = "eval_stats", .keptTrain = TRUE, ...){
  
  dotargs <- list(...)
  
  if(!.keptTrain){
    if(!("results" %in% names(wfRes))) 
      fold.res <- t(sapply(wfRes, function(x) 
        do.call(eval.function, c(list(trues=x$results$trues, 
                      preds=x$results$preds), dotargs))))
    else fold.res <- t(do.call(eval.function, c(list(trues = wfRes$results$trues, 
                                     preds = wfRes$results$preds), dotargs))) 
  }else{
    if(!("results" %in% names(wfRes))) 
      fold.res <- t(sapply(wfRes, function(x) 
        do.call(eval.function, c(list(trues=x$results$trues, 
                      preds=x$results$preds, 
                      y_train=x$train_y), dotargs))))
    else fold.res <- t(do.call(eval.function, c(list(trues = wfRes$results$trues, 
                                     preds = wfRes$results$preds, 
                                     y_train = wfRes$train_y), dotargs))) 
  }
  
  fold.res
}


#' Evalute the results of an internal predictive workflow
#' 
#' Calculate evaluation metrics from the raw results of a workflow
#' @param wfRes a data frame (or list of data frames) containing the results of
#' a predictive workflow with columns \code{trues} and \code{preds} containing
#' the real and predicted values, respectively
#' @param eval.function the function to be used to calculate error metrics from \code{wfRes}
#' @param y_train a vector of the whole training values
#' @param ... parameters to pass to \code{eval.function}
#'
#' @return The results (or a list of results) of \code{eval.function} applied to 
#' the data frame (or list of data frames) in \code{wfRes}
#' 
int_util_evaluate <- function(wfRes, eval.function = "eval_stats", y_train, ...){
    
  dotargs <- list(...)
  
    if(!("results" %in% names(wfRes))) 
      fold.res <- t(sapply(wfRes, function(x) 
        do.call(eval.function, c(list(trues = x$results$trues, 
                      preds = x$results$preds, 
                      y_train = y_train), dotargs))))
    else fold.res <- t( do.call(eval.function, c(list(trues = wfRes$results$trues, 
                                     preds = wfRes$results$preds, 
                                     y_train = y_train), dotargs)) ) 

  fold.res
}


#' Estimate error using a chosen method
#'
#' @param data a data frame
#' @param form a formula for learning
#' @param estimator the name of an error estimator function
#' @param est.pars a named list of arguments to feed to \code{estimator}
#' @param workflow the name of the workflow to use for making predictions
#' @param wf.pars a named list of arguments to feed to \code{workflow}
#' @param evaluator the name of the function to use to calculate evaluation results
#' @param eval.pars a named list of arguments to feed to \code{evaluator}
#' @param seed a seed to set before performing estimates
#'
#' @return The results of \code{evaluator} after applying \code{estimator} to the
#' learning task
#' 
#' @export
estimates <- function(data, form, estimator="kf_xval",
                      est.pars = list(nfolds=10, 
                                      fold.alloc.proc="Trand_SPrand"), 
                      workflow = "simple_workflow", wf.pars=NULL, 
                      evaluator = "evaluate", eval.pars=NULL,
                      seed=1234){
  
  if(!is.null(seed)) set.seed(seed)
  
  res <- do.call(estimator, c(list(data=data, form=form, 
                                   FUN=get(workflow, mode="function")), 
                              est.pars, wf.pars))
  
  est.res <- do.call(evaluator, c(list(wfRes=res), eval.pars))
  
  list(evalRes = est.res, rawRes = res)
}


############################
#### ESTIMATION METHODS
############################


#' Prequential evaluation
#'
#' Performs an evaluation procedure where training and test sets can 
#' be allocated in different ways, while always respecting the ordering 
#' provided by time (models are trained in the past and tested in the
#' relative future).
#' @param data full dataset
#' @param FUN function with arguments
#' \itemize{
#' \item \code{train} training set
#' \item \code{test} testing set
#' \item \code{time} column name of time-stamps
#' \item \code{site_id} column name of location identifiers
#' \item \code{form} a formula for model learning
#' \item \code{...} other arguments
#' }
#' @param form a formula for model learning
#' @param time column name of time-stamp in \code{data}. 
#' Default is "time"
#' @param site_id column name of location identifier in \code{data}. 
#' Default is "site_id"
#' @param .keepTrain if TRUE (default), instead of the results of 
#' \code{FUN} being directly returned, a list is created with both
#' the results and a \code{data.frame} with the time and site identifiers
#' of the observations used in the training step.
#' @param ... other arguments to FUN
#' @param nfolds number of folds for the data set to be separated into.
#' If you would like to set the number of time and space folds separately, 
#' \code{nfolds} should be set to \code{NULL} and \code{t.nfolds} and
#' \code{sp.nfolds} should be fed as a list to \code{alloc.pars}
#' (only available when using \code{fold.alloc.proc} \code{Tblock_SPrand})
#' @param window type of blocked-time window ordering considered. 
#' Should be one of
#' \itemize{
#' \item \code{growing} - for each time block being tested, all previous
#'  time blocks are used for training
#' \item \code{sliding} - for each time block being tested, the immediately
#'  previous time blocks are used for training
#' }
#' @param fold.alloc.proc name of fold allocation function. Should be one of
#' \itemize{
#' \item \code{Tblock_SPall} - each fold includes a block of contiguous time
#' for all locations
#' \item \code{Tblock_SPchecker} - each fold includes a block of contiguous time
#' for a systematically assigned (checkered) part of space
#' \item \code{Tblock_SPcontig} - each fold includes a block of contiguous time
#' for a block of spatially contiguous locations
#' \item \code{Tblock_SPrand} -  each fold includes a block of contiguous time
#' for a randomly assigned part of space
#' }
#' @param alloc.pars parameters to pass onto \code{fold.alloc.proc}
#' @param removeSP argument that determines whether spatio-temporal blocks
#' including the space being used for testing should be removed from the training set.
#' Default is FALSE, meaning the information is not removed
#' @param init_fold first temporal fold to use for testing. Default is 2.
#' @param .parallel Boolean indicating whether each block should be run in parallel
#' @param .verbose Boolean indicating whether updates on progress should be printed
#' @inherit kf_xval return
#' 
#' 
#' @import dplyr
#' @import foreach
#' 
#' @export
prequential_eval <- function(data, nfolds, FUN, form,
                             window = "growing", 
                             fold.alloc.proc="Tblock_SPall", alloc.pars=NULL, 
                             removeSP = FALSE, 
                             init_fold = 2,
                             time="time", site_id="site",
                             .keepTrain = TRUE, .parallel=TRUE, 
                             .verbose = ifelse(.parallel, FALSE, TRUE),
                             ...){
  requireNamespace("foreach", quietly=TRUE)
  
  assertthat::assert_that(is.data.frame(data),
                          time %in% colnames(data),
                          site_id %in% colnames(data),
                          window %in% c("growing", "sliding"),
                          fold.alloc.proc %in% c("Tblock_SPall",
                                                 "Tblock_SPchecker",
                                                 "Tblock_SPcontig",
                                                 "Tblock_SPrand"),
                          (is.null(alloc.pars) | is.list(alloc.pars)),
                          removeSP %in% c(TRUE, FALSE),
                          .keepTrain %in% c(TRUE, FALSE),
                          .parallel %in% c(TRUE, FALSE),
                          init_fold>=2, init_fold <= nfolds)
  
  `%myfun%` <- ifelse(.parallel, `%dopar%`, `%do%`)
  
  # call function that automatically allocates rows to folds
  fold_alloc <- do.call(fold.alloc.proc, c(list(data_coords = data[, c(time, site_id)], 
                                                nfolds = nfolds,
                                                time = time, 
                                                site_id = site_id), alloc.pars))
  data <- data[fold_alloc$idxs, ]
  folds <- fold_alloc$f
  
  assertthat::assert_that(is.vector(folds),
                          if(!is.null(nfolds)) length(unique(folds)) == nfolds)
  
  if(fold.alloc.proc != "Tblock_SPall")
    sep_folds <- apply(stringr::str_split_fixed(folds,"_", 2),
                       2, as.numeric)
  else
    sep_folds <- cbind(folds, 1)
  
  test_fold_inds <- which(sep_folds[,1] >= init_fold)
  test_folds <- sort(unique(folds[test_fold_inds]))
  
  if(.verbose) cat(paste0("Out of ", max(test_folds), ":"))
  
  pre.res <- foreach::foreach(f=test_folds) %myfun% {
    
    if(.verbose) cat(paste0(" ", f, "..."))
    
    if(fold.alloc.proc != "Tblock_SPall"){
      fs <- as.numeric(stringr::str_split_fixed(f,"_", 2))
      t.f <- fs[1]
      sp.f <- fs[2]
    }else{
      t.f <- as.numeric(f)
    }
    
    # each fold is used as test set once
    ts.inds <- which(folds == f)
    
    if(window=="growing"){
      if(!removeSP | fold.alloc.proc == "Tblock_SPall") 
        tr.inds <- which(sep_folds[,1] < t.f)
      else 
        tr.inds <- which(sep_folds[,1] < t.f &
                           sep_folds[,2] != sp.f)
    }else{
      if(!removeSP | fold.alloc.proc == "Tblock_SPall") 
        tr.inds <- which(sep_folds[,1] == (t.f-1))
      else 
        tr.inds <- which(sep_folds[,1] == (t.f-1) &
                           sep_folds[,2] != sp.f)
    }
    
    # split data into train and test sets
    train <- data[ tr.inds, ]
    test  <- data[ ts.inds, ]
    
    # call workflow returning results
    res <- FUN(form=form, train=train, test=test, time=time, site_id=site_id, ...)  
    if(.keepTrain) res <- c(res, list(train_y = train[ , as.character(form)[[2]]],
                                      train_idxs = fold_alloc$idxs[tr.inds],
                                      test_idxs = fold_alloc$idxs[ts.inds]))
    
    res
  }
  if(.verbose) cat("\n")
  
  pre.res
}

#' Cross-validation
#'
#' Performs a cross-validation experiment where folds can be 
#' allocated in different ways considering time and/or space
#' Performs a cross-validation experiment where folds can be 
#' allocated in different ways considering time and/or space
#' @param data full dataset
#' @param FUN function with arguments
#' \itemize{
#' \item \code{train} training set
#' \item \code{test} testing set
#' \item \code{time} column name of time-stamps
#' \item \code{site_id} column name of location identifiers
#' \item \code{form} a formula for model learning
#' \item \code{...} other arguments
#' }
#' @param form a formula for model learning
#' @param nfolds number of folds for the data set to be separated into. \cr
#' If you would like to set the number of time and space folds separately, 
#' \code{nfolds} should be set to \code{NULL} and \code{t.nfolds} and
#' \code{sp.nfolds} should be fed as a list to \code{alloc.pars}
#' (only available when using \code{fold.alloc.proc} set to 
#' \code{Tblock_SPchecker}, \code{Tblock_SPcontig} or \code{Tblock_SPrand}).
#' @param fold.alloc.proc name of fold allocation function. Should be one of
#' \itemize{
#' \item \code{Trand_SPrand} -- each fold contains completely random observations.
#' The default
#' \item \code{Tall_SPcontig} - each fold includes all time and a 
#' contiguous block of space
#' \item \code{Tall_SPrand} - each fold includes all time and random 
#' locations in space
#' \item \code{Tall_SPchecker} - each fold includes all time and a 
#' set of systematically assigned (checkered) part of space
#' \item \code{Tblock_SPall} - each fold includes a block of contiguous time
#' for all locations
#' \item \code{Trand_SPall} - each fold includes random time-snapshots of
#' of all locations
#' \item \code{Tblock_SPchecker} - each fold includes a block of contiguous time
#' for a systematically assigned (checkered) part of space
#' \item \code{Tblock_SPcontig} - each fold includes a block of contiguous time
#' for a block of spatially contiguous locations
#' \item \code{Tblock_SPrand} -  each fold includes a block of contiguous time
#' for a randomly assigned part of space
#' }
#' @param alloc.pars parameters to pass onto \code{fold.alloc.proc}
#' @param .parallel Boolean indicating whether each fold should be run in parallel
#' @inherit t_oos return
#' @param .verbose Boolean indicating whether updates on progress should be printed
#' @param time column name of time-stamp in \code{data}.  
#' @param site_id column name of location identifier in \code{data}. 
#' @param .keepTrain if TRUE (default), instead of the results of 
#' \code{FUN} being directly returned, a list is created with both
#' the results and a \code{data.frame} with the time and site identifiers
#' of the observations used in the training step.
#' @param ... other arguments to FUN
#' @return The results of \code{FUN}. Usually, a data.frame
#' with location identifier \code{site_id}, time-stamp \code{time},
#' true values \code{trues} and the workflow's predictions \code{preds}.
#' 
#' @import dplyr
#' @import foreach
#' 
#' @export
kf_xval <- function(data, nfolds, FUN, form,
                    fold.alloc.proc="Trand_SPrand", alloc.pars=NULL,
                    time="time", site_id="site",
                    .keepTrain=TRUE, .parallel=TRUE, 
                    .verbose = ifelse(.parallel, FALSE, TRUE), ...){
  
  requireNamespace("foreach", quietly=TRUE)
  
  assertthat::assert_that(is.data.frame(data),
                          time %in% colnames(data),
                          site_id %in% colnames(data),
                          fold.alloc.proc %in% c("Trand_SPrand",
                                                 "Tall_SPcontig",
                                                 "Tall_SPrand",
                                                 "Tall_SPchecker",
                                                 "Tblock_SPall",
                                                 "Trand_SPall",
                                                 "Tblock_SPchecker",
                                                 "Tblock_SPcontig",
                                                 "Tblock_SPrand"),
                          (is.null(alloc.pars) | is.list(alloc.pars)),
                          .keepTrain %in% c(TRUE, FALSE),
                          .parallel %in% c(TRUE, FALSE))
  
  # call function that automatically allocates rows to folds
  fold_alloc <- do.call(fold.alloc.proc, c(list(data_coords = data[, c(time, site_id)], 
                                                nfolds = nfolds,
                                                time = time, 
                                                site_id = site_id), alloc.pars))
  data <- data[fold_alloc$idxs, ]
  folds <- fold_alloc$f
  
  assertthat::assert_that(is.vector(folds),
                          if(!is.null(nfolds)) length(unique(folds)) == nfolds)
  
  `%myfun%` <- ifelse(.parallel, `%dopar%`, `%do%`)
  
  if(.verbose) cat(paste0("Out of ", max(folds), ":"))
  
  cv.res <- foreach::foreach(i=unique(folds)) %myfun% {
    
    if(.verbose) cat(paste0(" ", i, "..."))
    
    # each fold is used as test set once
    ts.inds <- which(folds == i)
    
    # split data into train and test sets
    train <- data[-ts.inds, ]
    test  <- data[ ts.inds, ]
    
    # call workflow returning results
    res <- FUN(form=form, train=train, test=test, time=time, site_id=site_id, ...)  
    if(.keepTrain) res <- c(res, list(train_y = train[ , as.character(form)[[2]]],
                                      train_idxs = fold_alloc$idxs[-ts.inds],
                                      test_idxs = fold_alloc$idxs[ts.inds]))
    res
  }
  
  if(.verbose) cat("\n")
  
  cv.res
}



############################
#### FOLD ALLOCATION METHODS
############################

#' Temporally blocked CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item blocked time
#' \item all locations
#' }
#' @param data_coords dataset with column time
#' @param nfolds number of folds
#' @param time column name of time-stamp in \code{data}. 
#' Default is "time"
#' @param site_id column name of location identifier in \code{data}. 
#' Default is "site"
#' @return A list with slots:
#' \itemize{
#' \item \code{idxs}, indices of data set as it should be re-ordered
#' \item \code{f}, a vector with the fold numbers 
#' (from \code{1} to \code{nfolds}) of each row in \code{data}
#' }
Tblock_SPall <- function(data_coords, nfolds, time="time", site_id="site"){
  
  # shuffle all rownames so that space will be shuffled
  shuffled.idxs <- shuffle(1:nrow(data_coords))
  
  # block time
  # sort time-stamps
  time.ids <- sort(unique(data_coords[[time]]))
  # divide time-stamps into folds
  t.folds <- cv_folds(time.ids, nfolds)
  names(t.folds) <- paste0("time_", time.ids)
  
  # assign rows to fold
  f <- t.folds[paste0("time_", data_coords[shuffled.idxs, time])]
  
  list(idxs = shuffled.idxs, f=f)
}
