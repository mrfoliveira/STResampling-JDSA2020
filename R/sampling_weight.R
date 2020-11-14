#' Calculate utility-based relevance
#' 
#' Calculate relevance of values given a parametrization
#' of the relevance function. 
#' Most relevant: phi -> 1; less relevant: phi -> 0.
#' @param y vector of values to calculate relevance of
#' @param phi.control list of parameters as returned
#' by function \code{UBL::phi.control}
#' @seealso \code{\link[UBL]{phi}}, \code{\link[UBL]{phi.control}}
get_phi <- function(y, phi.control){
  # require(UBL)
  UBL::phi(y, phi.control)
}

#' Calculate temporally-biased re-sampling weights without phi look-up
#'
#' Calculate weights for re-sampling with a temporal bias.
#' For each station, obsrvations that have been preceded by
#' extreme values the most recently, have higher weights
#' (meaning that observations where extreme values were recently
#' measured are more likely to be selected).
#' Most recent extreme observation in the past at station L: w -> 1; 
#' no rare observations in the past at station L: w -> 0.
#' 
#' @param times a vector of time-stamps
#' @param phi a vector of the relevance values of 
#' \code{df}'s target variable
#' @param rel.thr a relevance threshold above which an 
#' observation is considered relevant
#' 
#' @return A vector of temporally-biased re-sampling weights, scaled
#' to fit within range [0,1].
#' @author Mariana Oliveira
get_time_wts <- function(times, phi, rel.thr){
  # check types
  assertthat::assert_that(lubridate::is.Date(times) | lubridate::is.POSIXct(times) | 
                            lubridate::is.POSIXlt(times) | lubridate::is.POSIXt(times), 
                          msg = "times must be of type Date or POSIX")
 
  # scale time so most recent = 1
  time_wts <- as.numeric( lubridate::seconds( lubridate::interval(times, min(times)))) / 
    as.numeric( lubridate::seconds( lubridate::interval(max(times), min(times)) ))
 
  time_wts <- norm_scale(time_wts)
  
  time_wts
}


#' Calculate temporally-biased re-sampling weights with phi look-up
#'
#' Calculate weights for re-sampling with a temporal bias.
#' As time goes by past an extreme value, the weight tends to 1.
#' While no extreme value appears, the weight stays at 0 
#' (meaning they are less likely to be kept).
#' More distant to previous rare value: w -> 1; 
#' very close to recent rare value / no rare values so far: w -> 0.
#' 
#' @param df a data frame
#' @param phi a vector of the relevance values of 
#' \code{df}'s target variable
#' @param rel.thr a relevance threshold above which an 
#' observation is considered relevant
#' @param time the column name of the time-stamp
#' @param site_id the name of the column containing location IDs
#' 
#' @return A vector of temporally-biased re-sampling weights, scaled
#' to fit within range [0,1].
#' @author Mariana Oliveira
get_phi_time_wts <- function(df, phi, rel.thr, time, site_id){
  #require(dplyr, quietly=TRUE)
  
  times <- df[[time]]
  # check types
  assertthat::assert_that(lubridate::is.Date(times) | lubridate::is.POSIXct(times) | 
                            lubridate::is.POSIXlt(times) | lubridate::is.POSIXt(times), 
                          msg = "times must be of type Date or POSIX")
  
  new_df <- data.frame(time = df[[time]], station = df[[site_id]], phi = phi)
  time_diffs <- new_df %>% 
    dplyr::arrange(.data$time) %>%
    # see which cases are relevant
    dplyr::mutate(relev_flag = (.data$phi>=rel.thr)) %>%
    # calculate distances for each station separately
    dplyr::group_by(.data$station) %>% 
    dplyr::mutate(# new regime finds where regime changes from normal to extreme or vice-versa
                  # (start of series also considered "regime change" for this purpose)
                  new_regime = c(TRUE, (.data$relev_flag!=lag(.data$relev_flag))[-1]),
                  last_obs = .data$time==max(.data$time),
                  any_relev = any(.data$relev_flag),
                  # weight values at transition points
                  # (first relevant point, w = 1; first normal point, w=0)
                  val = ifelse(.data$new_regime & .data$relev_flag, 1,
                                  ifelse(.data$new_regime & !.data$relev_flag, 0 , NA)),
                  # if there is at least one rare value, then force most recent 
                  # observation to be considered "rare"; otherwise, keep it at 0
                  val = ifelse(.data$any_relev & .data$last_obs & !.data$new_regime, 1,
                               ifelse(!.data$any_relev & .data$last_obs, 0, .data$val)),
                  # interpolate between known values
                  interp = stats::approx(x=.data$time, y=.data$val, xout=.data$time)$y)

  # join to ensure correct row order
  new_df <- left_join(new_df, time_diffs, by = c("time", "station"))
  
  time_wts <- as.numeric(new_df[["interp"]])
  
  time_wts
}


#' Calculate spatially-biased re-sampling weights without phi look-up
#'
#' Calculate weights for re-sampling with a spatial bias without phi look-up.
#' Observations have a distance that tends to 1 as 
#' they are farther away from the spatially closest station
#' (meaning they are more likely to be kept).
#' Farthest away from other stations: d -> 1.
#' 
#' @param df a data frame
#' @param time the column name of the time-stamp
#' @param sites_sf An sf obejct containing station and IDs and 
#' geometry points of the locations. As an alternative, provide
#' \code{lon}, \code{lat}, and \code{crs}
#' @inheritParams df2site_sf
#'
#' @return A vector of spatially-biased re-sampling weights, scaled
#' to fit within range [0,1].
get_space_wts <- function(df, sites_sf=NULL, lon=NULL, lat=NULL, crs=NULL, site_id, time){
  
  if(is.null(sites_sf)){
    # get sites into right format
    assertthat::assert_that(!is.null(lon), !is.null(lat), !is.null(crs), 
                            msg = "Please provide locations object of type sf or 
                            CRS code and names of longitude and latitude columns")
    sites_sf <- df2site_sf(df, site_id, lon, lat, crs)
  }else{
    # check that stations in the data set match with stations in sites_sf
    assertthat::assert_that(all(df[[site_id]] %in% sites_sf[[site_id]]), 
                              msg = "please provide a locations object of type sf 
                            whose location names match with the names in df")  
  } 
  
  
  # create distance matrix
  dists <- get_spatial_dist_mat(sites_sf, site_id)
  
  # create data frame with minimum distances
  min_dists <- data.frame(station=unique(df[[site_id]]), min_dist=NA)
  for(i in 1:nrow(min_dists)){
    s <- min_dists[i,"station"]
    # check row for site s
    row <- which(rownames(dists)==paste0("SITE_",s))
    # check columns of sites that were relevant at this time slice (except itself)
    rmcol <- which(colnames(dists)==paste0("SITE_",s))  
    # minimum distance to other stations
    min_dists[i, "min_dist"] <- min(dists[row, -rmcol])
  }
  # rename rows for easy retrieval
  rownames(min_dists) <- paste0("SITE_", min_dists[,"station"])
  # normalize distance values so that max distance = 1
  min_dists[,"min_dist"] <- norm_scale(min_dists[,"min_dist"])
  
  # get space weights in the right order
  space_wts <- min_dists[paste0("SITE_", df[[site_id]]),"min_dist"]
  
  space_wts
}


#' Calculate spatially-biased re-sampling weights with phi look-up
#'
#' Calculate weights for re-sampling with a spatial bias with phi look-up.
#' Observations have a distance that tends to 1 as 
#' they are farther away from the closest relevant case (besides itself)
#' at time slice \code{t} (meaning they are more likely to be kept).
#' Farthest away from relevant cases at time slice t: d -> 1.
#' 
#' @param df a data frame
#' @param phi a vector of the relevance values of 
#' \code{df}'s target variable
#' @param rel.thr a relevance threshold above which an 
#' observation is considered relevant
#' @param time the column name of the time-stamp
#' @param sites_sf An sf obejct containing station and IDs and 
#' geometry points of the locations. As an alternative, provide
#' \code{lon}, \code{lat}, and \code{crs}
#' @inheritParams df2site_sf
#'
#' @return A vector of spatially-biased re-sampling weights, scaled
#' to fit within range [0,1].
get_phi_space_wts <- function(df, phi, rel.thr, sites_sf=NULL,
                          lon=NULL, lat=NULL, crs=NULL, site_id, time){
  
  if(is.null(sites_sf)){
    # get sites into right format
    assertthat::assert_that(!is.null(lon), !is.null(lat), !is.null(crs), 
                            msg = "Please provide locations object of type sf or 
                            CRS code and names of longitude and latitude columns")
    sites_sf <- df2site_sf(df, site_id, lon, lat, crs)
  }else{
    # check that stations in the data set match with stations in sites_sf
    assertthat::assert_that(all(df[[site_id]] %in% sites_sf[[site_id]]), 
                            msg = "please provide a locations object of type sf 
                            whose location names match with the names in df")  
  } 
  
  # create distance matrix
  dists <- get_spatial_dist_mat(sites_sf, site_id)
  max_dist <- max(dists)
  
  timz <- df[[time]]
  space_wts <- vector(mode="numeric", length=nrow(df))
  space_wts <- rep(NA, length(space_wts))
  for(i in 1:length(unique(timz))){
    # get time slice
    t <- unique(timz)[i]
    inds_t <- which(df[[time]]==t)
    
    # get indices of relevant cases at time slice t
    relev_inds <- inds_t[which(phi[inds_t] >= rel.thr)]
    
    # get indices of normal cases
    norm_inds <- setdiff(inds_t, relev_inds)
    if(!length(relev_inds)){
      # if there are no relevant cases, all have max distance 
      # (will be normalized to d=1)
      space_wts[inds_t] <- max_dist # 1
    }else{
      # otherwise, for each case 
      # find minimum distance to relevant case (at time slice t)
      relev_sites <- df[relev_inds, site_id]
      for(i in inds_t){
        s <- df[i, site_id]
        #if(length(setdiff(relev_sites, s))==0) browser()
        
        # if i is the only relevant case, it has maximum distance to other relevant cases
        if((length(unique(relev_sites))==1) && (s %in% relev_sites)){
          d <- max_dist # 1 
        # get minimum distance (to a relevant case)
        }else{
          # check row for site s
          row <- which(rownames(dists)==paste0("SITE_",s))
          # check columns of sites that were relevant at this time slice (except itself)
          cols <- which(colnames(dists) %in% paste0("SITE_", setdiff(relev_sites, s)))
          
          d <- min(dists[row, cols])
        } 
        
        # this is the raw space weight
        space_wts[i] <- d
      }
      
    }
    
    if(t==timz[1]) assertthat::assert_that(all(df[which(!is.na(space_wts)),time]==t))
  }
  
  # normalize weights to scale [0,1]
  space_wts <- norm_scale(space_wts)
  
  space_wts
}

#' Get spatio-temporal re-sampling weights
#'
#' A function that calculates different weights for
#' re-sampling that is temporally and/or spatially biased.
#' 
#' @details \code{phi} gives the target variable's relevance 
#' (higher relevance: phi -> 1; lower relevance: phi -> 0);
#' \code{time_wts} gives the observation's temporally biased
#' re-sampling weight (most recent observations: w -> 1; 
#' oldest: w -> 0.); \code{space_wts} gives the observation's
#' spatially biased re-sampling weight (farthest away from other 
#' relevant cases at time slice: d -> 1.).
#' High \code{time_wts} or \code{space_wts} means the observation is
#' more likely to be kept.
#' 
#' @param form a formula describing the learning task
#' @param df a data frame
#' @param alpha weighting parameter for temporal and spatial
#' re-sampling probabilities. Default 0.5
#' @param ptype type of spatio-temporal bias to calculate
#' @inheritParams get_phi
#' @inheritParams get_phi_space_wts
#' 
#' @return a data.frame with relevance \code{phi},
#' temporally biased weights \code{time_wts},
#' and spatially biased weights \code{space_wts} for
#' each row in \code{df}.
#'  
#' @seealso \code{\link{get_phi}}, \code{\link{get_time_wts}},
#'  \code{\link{get_phi_space_wts}}.
#'  
#' @export
sample_wts <- function(form, df, phi.control,
                       site_id, time, rel.thr, 
                       ptype, alpha,
                       sites_sf = NULL, lon=NULL, lat=NULL, crs = NULL){
  
  #----- ARGUMENT CHECK -------#
  
  assertthat::assert_that(xor(is.null(sites_sf), (is.null(lon) & is.null(lat) & is.null(crs))),
                          msg="Please provide either sites_sf or lon+lat+crs.")
  assertthat::assert_that(alpha>=0, alpha<=1, msg = "alpha must be between 0 and 1")
  assertthat::assert_that(time %in% colnames(df), site_id %in% colnames(df),
                          msg = "variables 'time' and 'site_id' must exist in df.")
  # check that there are no NAs in time and space tags
  assertthat::assert_that(!any(is.na(df[[time]])), 
                          !any(is.na(df[[site_id]])),
                          msg = "variables 'time' and 'site_id' cannot contain any NAs")
  assertthat::assert_that(ptype %in% c("orig", "addPhiLook", "addNoPhi"),
                          msg = "Please choose one of the available types of bias weight.")
  
  
  # check that there are no NAs in target
  y <- stats::model.response(stats::model.frame(form, df, na.action = NULL))
  assertthat::assert_that(!any(is.na(y)), 
              msg = "target variable must not contain any NAs")
  
  # check that either sites_sf or lon/lat are provided
  if(is.null(sites_sf)){
    assertthat::assert_that(!is.null(lon), !is.null(lat), !is.null(crs), 
                msg = "please provide locations object of type sf or 
                CRS code and names of longitude and latitude columns")
    assertthat::assert_that(!any(is.na(df[[lat]])), !any(is.na(df[[lon]])), 
                msg = "variables 'lat' and 'lon' cannot contain any NAs")
  }
  # check that stations in the data set match with stations in sites_sf
  if(!is.null(sites_sf)){
    assertthat::assert_that(all(df[[site_id]] %in% sites_sf[[site_id]]), 
                            msg = "please provide a locations object of type sf 
                            whose location names match with the names in df")
  }
  
  #----------------------------#
  
  
  # RELEVANCE
  phi <- get_phi(y, phi.control)
  
  # TIME
  timz <- df[[time]]
  if(alpha == 0){
    # if alpha = 0, only spatial weights are needed
    time_wts <- rep(0, nrow(df))
  }else{
    if(ptype=="addPhiLook"){
      time_wts <- get_phi_time_wts(df = df, phi = phi, rel.thr = rel.thr, 
                                       time = time, site_id = site_id)
    } else{
      time_wts <- get_time_wts(times = timz, phi = phi, rel.thr = rel.thr)
    }
      
  }
    
  # SPACE
  if(alpha==1){
    # if alpha = 1, only temporal weights are needed
    space_wts <- rep(0, nrow(df))
  }else{
    if(ptype == "addNoPhi")
      space_wts <- get_space_wts(df = df, 
                                 sites_sf = sites_sf, lon = lon, lat = lat, 
                                 crs = crs, site_id = site_id, time = time) 
    else
      space_wts <- get_phi_space_wts(df = df, phi = phi, rel.thr = rel.thr, site_id = site_id,
                                         sites_sf = sites_sf, lon = lon, lat = lat, time = time, crs = crs)
    
  }
  
  assertthat::assert_that(length(y)==length(phi),
              length(phi)==length(time_wts),
              length(time_wts)==length(space_wts))
  
  stprob <- data.frame(phi = phi, time_wts = time_wts, space_wts = space_wts,
                       stprob = alpha*time_wts + (1-alpha)*space_wts)
  
  # fix colnames
  if(ptype=="orig"){
    colnames(stprob)[c(3,4)] <- c("space_phi_wts", "stprob_orig")
  }else if(ptype=="addNoPhi"){
    colnames(stprob)[4] <- "stprob_addNoPhi"
  }else{
    colnames(stprob)[2:4] <- c("time_phi_wts", "space_phi_wts", "stprob_addPhiLook")
  }
    
  # get uncalculated weights back to NA
  if(alpha == 1)
    stprob[,3] <- as.numeric(NA)
  else if(alpha == 0)
    stprob[,2] <- as.numeric(NA)
  
  stprob
}
