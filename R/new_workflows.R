#' A simple learning and prediction workflow
#' 
#' A simple learning and prediction workflow that may deal
#' with NAs and use re-sampling techniques to balance an
#' imbalanced regression problem.
#' 
#' @param train a data frame for training
#' @param test a data frame for testing
#' @param time the name of the column in \code{train} and
#' \code{test} containing time-stamps
#' @param site_id the name of the column in \code{train} and
#' \code{test} containing location IDs
#' @param form a formula describing the model to learn
#' @param model the name of the algorithm to use
#' @param resample.pars parameters to be passed to re-sample function.
#' Default is NULL.
#' @param handleNAs string indicating how to deal with NAs.
#' If "centralImputNAs", training observations with at least 80\%
#' of non-NA columns, will have their NAs substituted by the mean
#' value and testing observatiosn will have their NAs filled in with
#' mean value regardless. Default is NULL.
#' @param min_train a minimum number of observations that must be
#' left to train a model. If there are not enough observations, 
#' predictions will be \code{NA}. Default is 2.
#' @param nORp a maximum number or fraction of columns/rows with missing
#' values above which a row/column will be removed from train before 
#' learning the model. Only works if \code{handleNAs} was
#' set to centralImputNAs. Default is 0.2.
#' @param ... other parameters to feed to \code{model}
#' 
#' @return a data frame containing time-stamps, location IDs,
#' true values and predicted values
#' 
#' @export
simple_workflow <- function(train, test, form, model,
                            time, site_id,
                            resample.pars = NULL,
                            handleNAs="centralImputNAs", 
                            min_train=2, nORp = 0.2, ...){
  
  #----- ARGUMENT CHECK -------#
  
  dotargs <- list(...)
  
  assertthat::assert_that(min_train>=2, 
                          msg = "Cannot train model with less than 2 observations.")
  assertthat::assert_that(time %in% colnames(train),
                            time %in% colnames(test),
                            msg = "'time' not a column in data set")
  assertthat::assert_that(site_id %in% colnames(train),
                            site_id %in% colnames(test),
                            msg = "'site_id' not a column in data set")
  
  
  #----------------------------#
  
  # get true values
  trues <- responseValues(form, test)
  tgt <- as.character(form[[2]])
  
  # save original state
  stats <- list()
  stats <- c(stats, list(
    nrow_tr_orig = nrow(train),
    trainCols_orig = colnames(train),
    tgt_sum_orig = summary(train[[tgt]])))
  
  col.inds <- which(colnames(train) %in% c(time, site_id))
  # correct default mtry if model is ranger and there is no argument given
  if(model=="ranger" & !("mtry" %in% names(dotargs)) & is.numeric(trues))
    dotargs$mtry <- max(floor(ncol(train[,-col.inds])/3), 1)
  
  # pre-process NAs
  if ( !is.null(handleNAs) & (anyNA(train) | anyNA(test)) ){
    if(handleNAs=="centralImputNAs"){
      test <- centralImputTsNAs(train, test, nORp)
      train <- centralImputTrNAs(train, nORp)
    }else{
      stop("different handling of NAs is future work")
    }
    stats <- c(stats, list(
      nrow_tr_postNA = nrow(train),
      trainCols_postNA = colnames(train),
      tgt_sum_postNA = summary(train[[tgt]])))
  }
  
  succ <- TRUE
  if(!is.null(resample.pars)){
    assertthat::assert_that(xor(is.null(resample.pars$C.perc), is.null(resample.pars$ratio)),
                            msg = "Please provide only one of either C.perc or ratio.")
    assertthat::assert_that( "type" %in% names(resample.pars), "thr.rel" %in% names(resample.pars),
                             msg = "Please provide type of resampling and relevance threshold.")
    
    ratio <- resample.pars$ratio
    # if ratio is provided, translate it into C.perc
    if(!is.null(ratio)){
      # save arguments
      rel <- resample.pars$rel
      # check defaults
      if(is.null(rel)) rel <- "auto"
      if(all(unlist(rel) == "auto")){
        cf <- if(!("cf" %in% names(resample.pars))) 1.5 else resample.pars$cf
      } 
      
      # calculate relevance
      relev <- calculate_relev(form = form, dat = train, rel = rel, cf = cf)
      # save calculated relevance to avoid extra calculations
      resample.pars$rel <- relev$pc
      # re-calculate C.perc
      C.perc <- tryCatch({ ratio2cperc(ratio = resample.pars$ratio, 
                           y.relev = relev$y.relev,
                           rel.thr = resample.pars$thr.rel,
                           type = resample.pars$type) },
                         error = function(err){
                           if(grepl("must be higher than the current ratio", err$message)){
                             warning(paste("Warning:", err$message), call. = FALSE)
                             succ <- FALSE
                             return(list(C.perc = NA, orig_ratio = NA))
                           }else{
                             stop(err)
                           }
                         })
        
      stats <- c(stats, list(tgt_ratio = ratio, orig_ratio = C.perc$orig.ratio,
                             calc_cperc = C.perc$C.perc))
      
      resample.pars$ratio <- NULL
      resample.pars$C.perc <- C.perc$C.perc
    }
    
    if(!is.na(resample.pars$C.perc)){
      tryCatch({ train <- do.call("RandomResample", c(list(form = form, 
                                                           dat = train, 
                                                           time = time, 
                                                           site_id = site_id),
                                                      resample.pars)) },
               error = function(err){
                 if(grepl("relevance function", err$message)){ 
                   warning("Relevance function should be redefined. Skipping resampling.", call. = FALSE)
                   succ <- FALSE
                 }else if(grepl("No examples will be removed or added", err$message)){
                   warning("C.perc is not useable. Skipping resampling.", call. = FALSE)
                   succ <- FALSE
                 }else{
                   stop(err) 
                 }
               })
    }
    
    succ <- succ & (nrow(train) != stats$nrow_tr_orig)
  }
  
  # check if experiment already failed
  if(!succ | nrow(train) < min_train){
    if(!succ) warning("Resampling failed. Predictions will be NA.", call. = FALSE)
    else warning("nrow(train)<min_train", call. = FALSE)
    
    res <- list( results = data.frame( trues=trues, preds=as.numeric(NA)) )
  }else{
    # check if columns need to be removed from train
    tr.col.inds <- which(colnames(train) %in% c(time, site_id))
    ts.col.inds <- which(colnames(test) %in% c(time, site_id))
    # removing offending column indices from train
    if(length(tr.col.inds)) train <- train[,-tr.col.inds]
    # train model
    m <- do.call(model, c(list(form, train), dotargs))
    # make predictions
    if(model=="ranger"){
      preds <- stats::predict(m, test[,-ts.col.inds])$predictions
    } else{
      preds <- stats::predict(m, test[,-ts.col.inds])
      if (is.numeric(train[[tgt]]) && !is.null(dim(preds)))
        preds <- preds[, 1]
    } 
    # prepare result object
    res <- list( results = data.frame(trues=trues, preds=preds) )
  }
  
  # save call
  args <- as.list(match.call())
  if(is.list(args$resample.pars)){
    if(!is.null(args$resample.pars$sites_sf))
      args$resample.pars$sites_sf <- NULL
  }
  
  stats <- c(stats, list(
    nrow_tr_final = nrow(train),
    trainCols_final = colnames(train),
    tgt_sum_final = summary(train[[tgt]])))

  res <- c(res,
           list(stats = stats),
           list(call. = args))
  
  res
}


#' A learning and prediction workflow with internal validation
#' 
#' A learning and prediction workflow that may deal
#' with NAs and use internal validation to parametrize a
#' re-sampling technique to balance an imbalanced regression problem.
#' 
#' @inheritParams simple_workflow
#' @param internal.est character string identifying the internal estimator
#'  function to use
#' @param internal.est.pars named list of internal estimator
#' parameters (e.g., tr.perc or nfolds)
#' @param internal.evaluator character string indicating internal 
#' evaluation function
#' @param internal.eval.pars named list of parameters to feed to internal 
#' evaluation function
#' @param stat parameter indicating summary statistic that should be 
#' used to determine the best internal evaluation metric: 
#' "MED" (for median) or "MEAN" (for mean)
#' @param metrics vector of names of two metrics to be used to determine the best 
#' parametrization (the second metric is only used in case of ties)
#' @param metrics.max vector of Booleans indicating whether each metric in 
#' parameter metrics should be maximized (TRUE) or minimized (FALsE)
#' for best results
#' @param resample.grid a data.frame with columns indicating
#' resample.pars to test using internal.est. Any NA value in resample.grid
#' will have the argument set to NULL.
#' @param .intRes a Boolean indicating whether the evalRes object outputed by
#' internal validation should be returned. Defaults to TRUE
#' @param .full_intRes a Boolean indicating whether the full results
#' object for internal validation should be returned as well. Defaults to FALSE
#' @param .int_parallel a Boolean indicating whether rows in the grid search
#' should be tested in parallel
#' 
#' @return a data frame containing time-stamps, location IDs,
#' true values and predicted values
#' 
#' @export
internal_workflow <- function(train, test, form, model, 
                              time, site_id,
                              resample.grid, 
                              resample.pars = NULL,
                              internal.est = NULL, 
                              internal.est.pars = NULL,
                              internal.evaluator = "int_util_evaluate", 
                              internal.eval.pars = NULL,
                              metrics = c("F1.u", "rmse_phi"), 
                              metrics.max = c(TRUE, FALSE), 
                              stat="MED",
                              handleNAs = "centralImputNAs", 
                              min_train = 2, nORp = 0.2,
                              .int_parallel = FALSE,
                              .intRes = TRUE,
                              .full_intRes = FALSE, ...){
  
  #----- ARGUMENT CHECK -------#
  
  dotargs <- list(...)
  
  
  assertthat::assert_that(min_train>=2,
                          msg = "Cannot train model with less than 2 observations.")
  assertthat::assert_that(time %in% colnames(train),
                            time %in% colnames(test),
                            msg = "'time' not a column in data set")
  assertthat::assert_that(site_id %in% colnames(train),
                            site_id %in% colnames(test),
                            msg = "'site_id' not a column in data set")
  assertthat::assert_that(!is.null(resample.pars),
                          msg = "Please provide minimum required resampling parameters.")
  assertthat::assert_that( "type" %in% names(resample.pars), "thr.rel" %in% names(resample.pars),
                           msg = "Please provide type of resampling and relevance threshold.")
  assertthat::assert_that( !("rel" %in% colnames(resample.grid)), !("cf" %in% names(resample.grid)),
                           msg = "rel and cf currently unavailable in resample.grid. (Future work!)")
  
  if(!all(unlist(resample.pars$rel) == unlist(internal.eval.pars$rel)) |
     !all(resample.pars$cf == internal.eval.pars$cf))
    warning("Be aware that relevance function arguments are different for resampling and evaluation.")
  
  # check relevance arguments
  rs.rel <- resample.pars$rel
  ev.rel <- internal.eval.pars$rel
  # check rel arguments
  if (is.null(rs.rel) & !is.null(ev.rel)){
    warning("Setting resample.pars$rel to internal.eval.pars$rel.")
    rs.rel <- ev.rel
  }else if (!is.null(rs.rel) & is.null(ev.rel)){
    ev.rel <- rs.rel 
    warning("Setting internal.eval.pars$rel to resample.pars$rel.")
  }else if (is.null(rs.rel) & is.null(ev.rel)){
    rs.rel <- ev.rel <- "auto" 
  }
  
  rs.cf <- resample.pars$cf
  ev.cf <- internal.eval.pars$cf
  if (is.null(rs.cf) & !is.null(ev.cf) & all(unlist(rs.rel) == "auto")){
    warning("Setting resample.pars$cf to internal.eval.pars$cf.")
    rs.cf <- ev.cf
  }else if (!is.null(rs.cf) & is.null(ev.cf) & all(unlist(ev.rel) == "auto")){
    ev.cf <- rs.cf 
    warning("Setting internal.eval.pars$cf to resample.pars$cf.")
  }else if (is.null(rs.cf) & is.null(ev.cf)){
    if(all(unlist(rs.rel) == "auto"))
      rs.cf <- 1.5
    if(all(unlist(ev.rel) == "auto"))
      ev.cf <- 1.5
  }
  
  #----------------------------#
  
  trues <- responseValues(form, test)
  tgt <- as.character(form[[2]])
  col.inds <- which(colnames(train) %in% c(time, site_id))
  
  # calculate phi.control with whole training set so that both resampling
  # and internal evaluation is done with the same relevance function
  if(all(unlist(rs.rel) == unlist(ev.rel)) & rs.cf == ev.cf){
    rs.rel <- ev.rel <- get_phi_control(y = train[,tgt], rel = rs.rel, cf = rs.cf)
  }else{
    rs.rel <- get_phi_control(y = train[,tgt], rel = rs.rel, cf = rs.cf)
    ev.rel <- get_phi_control(y = train[,tgt], rel = ev.rel, cf = ev.cf)
  }

  # saving the results of phi.control
  resample.pars$rel <- rs.rel
  int_resample.pars <- resample.pars
  internal.eval.pars$rel <- ev.rel
  
  assertthat::assert_that(is.list(resample.pars$rel),
                          msg = "Something went wrong when calculating relevance.")
  
  # the relevance function for internal evaluation should be calculated 
  # using the whole available training values when calculating internal evaluation metrics
  if(!("y_train" %in% names(internal.eval.pars)))
    internal.eval.pars$y_train <- train[, tgt]
  
  # correct default mtry if model is ranger and there is no argument given
  if(model=="ranger" & !("mtry" %in% names(dotargs)) & is.numeric(responseValues(form, test)))
    dotargs$mtry <- max(floor(ncol(train[,-col.inds])/3), 1)
  
  # prepare data frame for grid results
  metrics.df <- data.frame(matrix(ncol = length(metrics)*5, nrow = nrow(resample.grid)))
  colnames(metrics.df) <- c(paste0("MED", metrics), paste0("IQR", metrics),
                            paste0("MEAN", metrics), paste0("SD", metrics),
                            paste0("nonNAs", metrics))
  grid.res <- cbind(resample.grid, metrics.df)
  
  `%myfun%` <- ifelse(.int_parallel, `%dopar%`, `%do%`)
  # do internal validation 
  int.results <- foreach::foreach(i = 1:nrow(resample.grid)) %myfun% {
  # int.results <- list()
  # for(i in 1:nrow(resample.grid)){
     
    # update resample.pars to params to test internally
    for(j in 1:ncol(resample.grid)){
      if(!is.na(resample.grid[i,j])){
        # save argument
        int_resample.pars[[colnames(resample.grid)[j]]] <- resample.grid[i,j] }
      else{
        # argument is discarded
        int_resample.pars[[colnames(resample.grid)[j]]] <- NULL
      } 
    }
      
    # get results with these parameters
    int.res <- NULL
    int.res <- estimates(data = train, form = form, 
                         estimator = internal.est,
                         est.pars = internal.est.pars, 
                         workflow = "simple_workflow", 
                         wf.pars = c(list(model = model,
                                          handleNAs = handleNAs, 
                                          min_train = min_train, 
                                          nORp = nORp,
                                          resample.pars = int_resample.pars), 
                                     dotargs), 
                         evaluator = internal.evaluator, 
                         eval.pars = internal.eval.pars, 
                         seed = NULL)
    
    i_results <- list()
    # if re-sampling succeeded, save results
    if(.full_intRes){
      i_results <- c(int.res, list(int_call. = int.res$rawRes[[1]]$call.))
    }else{
      i_results <- list(evalRes = int.res$evalRes, 
                          int_call. = int.res$rawRes[[1]]$call.)
    } 
      
    grid.res[i, c((ncol(resample.grid) + 1):ncol(grid.res))] <- summarize_metrics(int.res$evalRes, metrics)
    
    i_results
    # int.results[[i]] <- i_results
  }
  
  if(all(is.na(grid.res[, paste0(stat, metrics)]))){
    warning("Internal evaluation failed. Predictions will be NA.")
    
    res <- list(results = data.frame(trues=trues, preds=as.numeric(NA)),
                resample_grid = grid.res,
                resample_chosen = resample.grid[integer(0),])
    
  }else{

    best.idx <- get_best(grid.res, metrics, metrics.max, stat)
    # update resample.pars to best results
    for(j in 1:ncol(resample.grid)){
      if(!is.na(resample.grid[best.idx,j])){
        # save argument
        resample.pars[[colnames(resample.grid)[j]]] <- resample.grid[best.idx,j] }
      else{
        # argument is discarded
        resample.pars[[colnames(resample.grid)[j]]] <- NULL
      } 
    }
    
    res <- do.call("simple_workflow", c(list(train = train, 
                                             test = test, form = form,
                                             model = model,
                                             time = time, site_id = site_id,
                                             handleNAs = handleNAs, 
                                             min_train = min_train, 
                                             nORp = nORp,
                                             resample.pars = resample.pars),
                                        dotargs))
    
    # save more info
    res <- c(res, list(resample_grid_results = grid.res,
                resample_chosen = resample.grid[best.idx,]))
    names(res)[which(names(res) == "call.")] <- "internal_call."
    res$internal_call.$train <- NULL
    res$internal_call.$test <- NULL
   }

  # save call
  args <- as.list(match.call())
  if(is.list(args$resample.pars)){
    if(!is.null(args$resample.pars$sites_sf))
      args$resample.pars$sites_sf <- NULL
  }
  
  # save call
  res <- c(res, list(call. = args))
  # save internal results
  if(.intRes) res <- c(res, list(internal_results = int.results))
  
  res
}
