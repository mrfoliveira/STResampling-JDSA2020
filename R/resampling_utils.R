#' Calculate C.perc from ratio
#'
#' @param ratio Value that defines the target ratio between the number of 
#' extreme values and the number of normal values. If more than one bump exists:
#' re-sampling will be done proportionally to each bump, so that the ratio between
#' the sum of all extreme values and all normal values in the final data set
#' is the target ratio.
#' @param y.relev A vector containing the relevance values of target
#' @param rel.thr A relevance threshold
#' @param type "under", "over", or "gauss"
#'
#' @return A list with C.perc, the re-sampling percentage that will result 
#' in the target ratio and orig_ratio indicating the original ratio between
#' extreme and normal values
#' 
#' @export
ratio2cperc <- function(ratio, y.relev, rel.thr, type){
  assertthat::assert_that(type %in% c("under", "over", "gauss"))
  assertthat::assert_that(is.numeric(ratio), length(ratio)==1,
                          msg="Please provide only one numeric extreme/normal ratio")
  assertthat::assert_that(ratio > 0,
                          msg = "Target extreme/normal ratio must be over 0")
  
  extr.counts <- length(which(y.relev >= rel.thr))
  norm.counts <- length(y.relev) - extr.counts
  
  cur.ratio <- extr.counts / norm.counts
  
  assertthat::assert_that(ratio > cur.ratio, 
                          msg = "The target ratio (extreme/normal) must be higher than the current ratio.")
  
  ttype <- ifelse(type=="under", "und", "ove")
  C.perc <- switch (ttype,
                    "und" = extr.counts / (norm.counts * ratio),
                    "ove" = (ratio*norm.counts - extr.counts) / extr.counts )
  C.perc <- round(C.perc, 5)
  
  if(type == "under"){
    if(round(C.perc*norm.counts) <1){
      warning("The target ratio is too close to the current ratio.", call. = FALSE)
      C.perc <- NA
    }
  }else{
    if(round(C.perc*extr.counts) <1 | (type == "gauss" & round(C.perc*extr.counts) <2)){
      warning("The target ratio is too close to the current ratio.", call. = FALSE)
      C.perc <- NA
    }
  }
  
  list(C.perc = C.perc, orig.ratio = cur.ratio)
}

#' Check and re-calculate C.perc if necessary
#'
#' @param C.perc Value or vector of values with re-sampling percentages. 
#' @param inds A list with entries obs.ind (a list with named vector cointaining
#' target values for each bump), und, and ove (containing the indices of 
#' obs.ind that correspond to bumps of extreme and normal values, respectively).
#' @param type "under", "over" or "gauss"
#'
#' @return A vector of re-sampling percentages for each bump that will be 
#' subject to re-sampling
check_cperc <- function(C.perc, inds, type){
  
  norm.counts <- sapply(inds$obs.ind[inds$und], length)
  extr.counts <- sapply(inds$obs.ind[inds$ove], length)
  
  ttype <- ifelse(type=="under", "und", "ove")
  if (is.vector(C.perc) && is.numeric(C.perc)) {
    if (length(inds[[ttype]]) > 1 & length(C.perc) == 1) {
      C.perc <- rep(C.perc[1], length(inds[[ttype]]))
    }
  } 
  
  assertthat::assert_that( length(inds[[ttype]]) == length(C.perc) ,
                           msg = "The number of re-sampling percentages must be one or correspond to the number of bumps below the threshold defined!")
  assertthat::assert_that(is.vector(C.perc), is.numeric(C.perc), 
                          msg = "C.perc must be a number or vectors of numbers")
  assertthat::assert_that(all(C.perc > 0),
                          msg = "C.perc must be over 0")
  if(type == "under"){
    assertthat::assert_that(all(C.perc < 1),
                            msg = "Under-sampling percentage must be below 1")
    assertthat::assert_that( any(round(C.perc*norm.counts) >0),
                             msg = "No examples will be removed or added. 
                             Please redefine your arguments.")
  }
  else if(type == "over")
    assertthat::assert_that( any(round(C.perc*extr.counts) >0),
                             msg = "No examples will be removed or added. 
                             Please redefine your arguments.")
  else
    assertthat::assert_that( any(round(C.perc*extr.counts) >1),
                             msg = "No examples will be removed or added. 
                             Please redefine your arguments.")
  
  C.perc
}

#' Calculate bias weights
#' 
#' @inheritParams sample_wts
#' @param dat A data frame with time and site_id columns.
#' @param thr.rel a relevance threshold above which an 
#' observation is considered relevant.
#' @param ptype Type of spatio-temporal bias to calculate.
#' Can include "orig", "addPhiLook", "addNoPhi"
#' @param ... Parameters to pass to 
#' \code{sample_wts}. Can include weight type (defaults to orig), 
#' alpha (defaults to 0.5), and location info (sites_sf or lat+lon+crs).
#'
#' @return Vector of bias weights
calculate_bias <- function(form, dat, phi.control, thr.rel, 
                           time, site_id, ptype, ...){
  
  # calculate sample weights
  stprobs <- sample_wts(form = form, df = dat, 
                        phi.control = phi.control,
                        rel.thr = thr.rel, 
                        time = time, site_id = site_id,
                        ptype = ptype, ...)
  stprobs <- stprobs[ , paste0("stprob_", ptype), drop=TRUE]
  names(stprobs) <- rownames(dat)
  
  stprobs
}

#' Select cases for adding / keeping in data set
#'
#' @param obs.ind Named vector of target values
#' @param C.perc Percentage of cases to select from obs.ind
#' @param repl  A Boolean value controlling whether replication is allowed
#' when re-sampling observations. Defaults to FALSE when under-sampling
#' and to TRUE when over-sampling or adding Gaussian noise.
#' @param prob Named vector of selection probability. Defaults to NULL.
#' @param epsilon Mimimum probability value to add to prob after normalization.
#'
#' @return Vector containing selected observation names 
select_cases <- function(obs.ind, C.perc, repl, prob=NULL, epsilon=1E-4){
  if(!is.null(prob)){
    # normalize just cases in this bump and add epsilon
    prob <- prob[match(names(obs.ind), names(prob))]
    prob <- norm_scale(prob) + epsilon 
  }
  
  # sample indices of examples to keep/add
  sel <- sample(x = names(obs.ind), 
                size = round(C.perc * length(obs.ind)), 
                replace = repl, 
                prob = prob)
  
  sel
}


#' Provide examples with added Gaussian noise
#'
#' @param exs Data frame with observations to perturb
#' @param bump.exs Data frame with observations to use as baseline sd
#' @param exclude Indices of columns to exclude from Gaussian noise 
#' (e.g., spatio-temporal coordinates)
#' @param pert Percentage of the sd of bump.exs that should be used as
#' sd for the gaussian noise added to exs
#'
#' @return Data frame where all numeric values have added Gaussian noise
#' (except for the columns in exclude)
gauss.exs <- function(exs, bump.exs, pert, exclude=integer(0)){
  # save column name order
  cnames <- colnames(exs)
  
  # save columns to exclude from noise
  if(length(exclude)){
    exs.excl <- exs[ , exclude]
    exs <- exs[ , -exclude]
    bump.exs <- bump.exs[ , -exclude]
  }
  
  for(j in 1:ncol(exs)){
    if(!is.numeric(bump.exs[ ,j])){
      warning("Non-numeric column: did not add Gaussian noise.")
    }else{
      # calculate sd based on the whole bump
      pert_sd <- pert * stats::sd(bump.exs[ , j], na.rm=TRUE)
      assertthat::assert_that( !is.na(pert), msg = "Calculating sd failed.")
      
      # perturb numeric columns
      exs[ ,j] <- exs[ ,j] + stats::rnorm(1, 0, pert_sd)
    }
  }
  
  # join with unperturbed columns in the original column order
  if(length(exclude)>0) exs <- cbind(exs.excl, exs)[, cnames]
  
  exs
}

#' Get phi control
#'
#' @param y the target variable
#' @inheritParams RandomResample
#'
#' @return phi.control object
get_phi_control <- function(y, rel, cf = 1.5){
 
  if (is.matrix(rel)) {
    pc <- uba::phi.control(y, method = "range", control.pts = rel)
  }
  else if (is.list(rel)) {
    pc <- rel
  }
  else if (rel == "auto") {
    pc <- uba::phi.control(y, method = "extremes", coef = cf)
  }
  else {
    stop("Argument rel should be 'auto', a list returned by phi.control or a control.pts matrix.")
  }
  pc
}

#' Calculate phi control
#'
#' @inheritParams RandomResample
#'
#' @return A list containing a vector s.y containing sorted target values,
#' a vector y.relev containing the relevance values of s.y, and phi.control
calculate_relev <- function(form, dat, rel, cf = 1.5){
  
  tgt <- which(names(dat) == as.character(form[[2]]))
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  # get arguments
  pc <- get_phi_control(y, rel, cf = 1.5)
  y.relev <- UBL::phi(s.y, pc)
  
  list(s.y=s.y, y.relev=y.relev, pc = pc)
}

#' Calculate bumps
#'
#' @param y.relev Vector containing relevance of sorted target values
#' @inheritParams RandomResample
#'
#' @return Vector containing the first index of each bump
calculate_bumps <- function(y.relev, thr.rel){
  bumps <- c()
  for (i in 1:(length(y.relev) - 1)) {
    if ((y.relev[i] >= thr.rel && y.relev[i + 1] < thr.rel) || 
        (y.relev[i] < thr.rel && y.relev[i + 1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  
  bumps
}

#' Calculate list of indices in each bump
#'
#' @param s.y Vector of sorted target values
#' @param bumps Vector containing first index of each bump
#'
#' @return List where each entry has a named vector of target values
#' belonging to a single bump
calculate_obs_ind <- function(s.y, bumps){
  
  nbump <- length(bumps) + 1
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]
  
  obs.ind
}

#' Calculate indices of relevance bumps
#'
#' @param s.y Vector of sorted target values
#' @param y.relev Vector containing relevance of sorted target values
#' @param pc Relevance phi.control
#' @inheritParams RandomResample
#'
#' @return A list with entries obs.ind (a list with named vector cointaining
#' target values for each bump), und, and ove (containing the indices of 
#' obs.ind that correspond to bumps of extreme and normal values, respectively)
get_inds <- function(s.y, y.relev, pc, thr.rel){
  
  # get relevance
  assertthat::assert_that(any(y.relev<1), any(y.relev>0),
                          msg="All points have the same (0 or 1) relevance. Please, redefine your relevance function.")
  assertthat::assert_that(any(y.relev >= thr.rel), any(y.relev < thr.rel),
                          msg = "There must be at least one normal and one extreme value in the data set. Please, redefine your relevance function.")
  
  # calculate bumps in relevance
  bumps <- calculate_bumps(y.relev, thr.rel)
  
  # get list with relevant indices
  obs.ind <- calculate_obs_ind(s.y, bumps)
  imp <- sapply(obs.ind, function(x) mean(UBL::phi(x, pc)))
  und <- which(imp < thr.rel)
  ove <- which(imp >= thr.rel)
  
  list(obs.ind = obs.ind, und = und, ove = ove)
}