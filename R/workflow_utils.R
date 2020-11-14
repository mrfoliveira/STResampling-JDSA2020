
#' Get index of best results
#'
#' @param grid.res data frame with parameters and metrics results
#' @param metrics vector of (up to 2) metrics that should be 
#' used to decide best internal results
#' @param metrics.max vector of Boolean value indicating whether metrics 
#' should be maximized for best results (in corresponding order to metrics)
#' @param stat Stat that should be used to decide on best result
#' (can be "MED" (median), "MEAN", "IQR" or "SD")
#'
#' @return the row index of the best results
get_best <- function(grid.res, metrics, metrics.max, stat){
  # find best results from internal validation
  if(metrics.max[1] == TRUE){
    best1 <- max(grid.res[,paste0(stat, metrics[1])], na.rm=TRUE) 
  }else{
    best1 <- min(grid.res[,paste0(stat, metrics[1])], na.rm=TRUE) 
  }
  best <- which(grid.res[,paste0(stat, metrics[1])] == best1)
  
  # handle ties
  if(length(best) >1 & length(metrics) >1){
    if(metrics.max[2] == TRUE)
      best2 <- max(grid.res[best, paste0(stat, metrics[2])], na.rm=TRUE)
    else
      best2 <- min(grid.res[best, paste0(stat, metrics[2])], na.rm=TRUE)
    
    best <- which(grid.res[ , paste0(stat, metrics[1])] == best1 & 
                    grid.res[ , paste0(stat, metrics[2])] == best2)
  }
  
  if(length(best) >1){
    warning("Could not resolve metrics tie. Using the first set of best parameters.")
    best <- best[1]
  }
  
  best
}

#' Summarize metrics
#'
#' Used for internal validation
#' 
#' @param int.res A list with results obtained by running estimates
#' @param metrics a list of metrics that should be summarized
#'
#' @return a named vector with median, IQR, mean, standard-deviation,
#' and number of non-NA values of each metric
summarize_metrics <- function(int.res, metrics){
  
  int.res <- int.res[, metrics, drop=F]
  
  medians <- apply(int.res, 2, stats::median, na.rm=TRUE)
  iqrs <- apply(int.res, 2, stats::IQR, na.rm=TRUE)
  means <- apply(int.res, 2, mean, na.rm=TRUE)
  sds <- apply(int.res, 2, stats::sd, na.rm=TRUE)
  nonNAs <- apply(int.res, 2, function(x) length(which(!is.na(x))))
  
  nms <- c(paste0("MED", metrics), paste0("IQR", metrics),
           paste0("MEAN", metrics), paste0("SD", metrics),
           paste0("nonNAs", metrics))
  
  res <- c(medians, iqrs, means, sds, nonNAs)
  names(res) <- nms
  res
}



#' Handling NAs in train and test
#' 
#' Discard columns/rows with too many NAs and 
#' then impute with central value.
#'
#' @param train training data set
#' @param test testing data set
#' @param nORp minimum percentage of NA values in a row/column for
#' that row/column to be discarded from training set
#'
#' @return list with an entry for the training and test sets
#' (in this order), both now with no NA values
#' 
#' @export
centralImputTsNAs <- function(train, test, nORp){
  
  if(anyNA(test)){
    # fill in test columns with central value of train column
    for (i in seq(ncol(test))) if (any(idx <- is.na(test[, i]))){
      test[idx, i] <- DMwR2::centralValue(train[, i])          
    } 
  }
  
  test
}


#' Handling NAs in train and test
#' 
#' Discard columns/rows with too many NAs and 
#' then impute with central value.
#'
#' @param train training data set
#' @param nORp minimum percentage of NA values in a row/column for
#' that row/column to be discarded from training set
#'
#' @return list with an entry for the training and test sets
#' (in this order), both now with no NA values
#' 
#' @export
centralImputTrNAs <- function(train, nORp){
  # discard columns with too many NAs if there would be still predictors left
  discCols <- which(sapply(train, 
                           function(y) length(which(is.na(y)))/length(y)) > nORp)
 
  if(length(discCols) > 0) {
    if( (ncol(train) - length(discCols) - 3) > 0 ){
      warning(paste("Dropped", length(discCols), "columns from train due to NAs. Keeping", ncol(train) - length(discCols),"variables."), call. = FALSE)
      train <- train[, -discCols]
    }else{
      warning(paste("Should have dropped", length(discCols), "columns from train due to NAs. Keeping all because too few would be left."), call. = FALSE)
    }
  }
  
  # discard rows with too many NAs in train
  suppressWarnings( idxs <- DMwR2::manyNAs(train, nORp = nORp) )
  if(length(idxs)){
    warning(paste0("Dropped ", length(idxs), " rows from train due to NAs. Keeping ", nrow(train) - length(idxs)," (", round(100*(nrow(train) - length(idxs))/nrow(train)),"%) rows."), call. = FALSE)
    train <- train[-idxs, ]
  } 
  
  # fill in empty value in train
  if(anyNA(train)) train <- DMwR2::centralImputation(train)
  
  train
}
