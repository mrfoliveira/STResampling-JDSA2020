#' Random re-sampling for imbalanced regression problems (with spatio-temporal bias)
#' 
#' @param form A formula describing the prediction problem.
#' @param dat A data frame containing the original imbalanced data set.
#' @param type Type of re-sampling to apply. Can be one of "under", "over", 
#' and "gauss", depending on whether the user wants to under-sample normal cases,
#' over-sample extreme cases or add Gaussian noise to replicated extreme cases.
#' @param rel The relevance function which can be automatically ("auto") 
#' determined (the default) or may be provided by the user through a matrix
#' with the interpolating points.
#' @param thr.rel A number indicating the relevance threshold below which a 
#' case is considered as belonging to the normal "class".
#' @param cf Parameter needed if rel = 'auto'. The default is 1.5.
#' @param C.perc Vector containing percentage values (or a single  value that 
#' will be used for all bumps. In under-sampling, C.perc of the size of each 
#' bump of normal values will be kept in the final data set. In the case of 
#' over-sampling and Gaussian noise, C.perc of the size of each bump of extreme
#' values will be added to the final data set. Bumps are ordered in ascending 
#' order of the target value.
#' @param repl A Boolean value controlling whether replication is allowed
#' when re-sampling observations. Defaults to FALSE when under-sampling
#' and to TRUE when over-sampling or adding Gaussian noise.
#' @param pert Standard deviation of gaussian noise as a percentage of
#' of each variable original standard deviation. Only necessary if type = "gauss"
#' @param bias Boolean indicating whether spatio-temporal bias should be
#' factored in while re-sampling
#' @param time Column name of the time-stamp (if available). 
#' Only necessary if bias = TRUE or type = "gauss"
#' @param site_id Column containing location IDs (if available).
#' Only necessary if bias = TRUE or type = "gauss"
#' @param ... Parameters to feed to \code{sample_wts} in case bias = TRUE.
#'
#' @references Paula Branco, Rita P. Ribeiro, Luis Torgo (2016)., 
#' UBL: an R Package for Utility-Based Learning, 
#' CoRR abs/1604.08079 [cs.MS], URL: http://arxiv.org/abs/1604.08079
#' 
#' @seealso \code{\link[UBL]{RandUnderRegress}}, \code{link{sample_wts}}.
#' 
#' @export
RandomResample <- function (form, dat, type, C.perc,
                            thr.rel, rel = "auto", cf=1.5,
                            repl = ifelse(type=="under", FALSE, TRUE), 
                            pert = 0.1, time = NULL, site_id = NULL,
                            bias = FALSE, ...) 
{

  #----- ARGUMENT CHECK -------#
  
  # check main arguments
  assertthat::assert_that(type %in% c("under", "over", "gauss"))
  
  # check arguments for bias
  if(bias)
    assertthat::assert_that(!is.null(time), !is.null(site_id),
                            msg = "Please make sure the required parameters for biased re-sampling are provided.") 
  
  # check arguments for gaussian
  if(type=="gauss"){
    # check columns that should be excluded from noise
    if(!is.null(time) | !is.null(site_id)){
      coords <- which(colnames(dat) %in% c(time, site_id))
      assertthat::assert_that(length(coords)>0, 
                              msg = "The coordinate names time and/or site_id provided do not correspond to columns in the data.")
      assertthat::assert_that(any(unlist(sapply(dat[, -coords], is.numeric))),
                              msg = "There are no numeric (non-coordinate) variables that can have gaussian noise addded to them.")
    }else{
      coords <- integer(0) 
    }
  }
  
  #------- RE-SAMPLING -------#
  
  # calculate relevance values
  relev <- calculate_relev(form = form, dat = dat, rel = rel, cf = cf)
    
  # calculate indices for each bump
  inds <- get_inds(s.y = relev$s.y, y.relev = relev$y.relev, 
                   pc = relev$pc, thr.rel = thr.rel)
  
  # make sure C.perc works for this data
  C.perc <- check_cperc(C.perc = C.perc, inds = inds, type = type)
  
  # calculate bias probabilities
  if(bias){   
    bias.params <- list(...)
    if("epsilon" %in% names(bias.params)){
      epsilon <- bias.params$epsilon
      bias.params$epsilon <- NULL
    }else{
      warning("Undefined epsilon value. Using 1E-4.")
      epsilon <- 1E-4
    } 
  
    bias.prob <- do.call("calculate_bias", c(list(form = form, dat = dat, 
                              phi.control = relev$pc,
                              thr.rel = thr.rel, 
                              time = time, site_id = site_id), bias.params))
  
  } else{
    bias.prob <- NULL
    epsilon <- NULL
  }
  
  # get type of bumps
  ttype <- ifelse(type=="under", "und", "ove")
  if(type=="under"){
    # keep all extreme cases
    newdata <- NULL
    for (j in 1:length(inds$ove))
      newdata <- rbind(newdata, dat[names(inds$obs.ind[[inds$ove[j]]]), ])
    
  }else{
    # check argument
    if(any(C.perc >1))
      assertthat::assert_that(repl,
                              msg = "When over-sampling percentage is above 1, repl should be set to TRUE")
    
    # keep all cases
    newdata <- dat 
  }
    
  # keep a percentage of normal cases (under)
  # or add a percentage of extreme cases (over)
  for (j in 1:length(inds[[ttype]])) {
    b <- inds[[ttype]][j]
    
    # get indices of cases to add/keep
    sel <- select_cases(obs.ind = inds$obs.ind[[b]],
                        C.perc = C.perc[[j]], 
                        repl = repl, prob = bias.prob, epsilon = epsilon)
    new.exs <- dat[sel, ] 
    
    if(length(sel) == 0){
      warning(paste("No examples selected for bump", j))
    }else{
      # add noise to examples
      if(type=="gauss"){
        if(length(sel) >1){
          # all examples in this bump
          bump.exs <- dat[names(inds$obs.ind[[b]]), ]
          # coordinates should be kept unperturbed
          new.exs <- gauss.exs(exs = new.exs, bump.exs = bump.exs, pert = pert, exclude = coords) 
        }else{
          warning(paste0("Only one example in bump ", j,". Noise was not added."))
        }
      }
    }
    
    # add examples to data
    newdata <- rbind(newdata, new.exs)
  }
  newdata
}