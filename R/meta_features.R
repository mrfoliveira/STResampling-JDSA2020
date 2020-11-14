#########

get_stats <- function(x, vname="", sum_stats=c("min", "median", "quantile", "mean", "max", "sd"),
                      probs = c(0.25, 0.75), na.rm=TRUE){
  x <- x[!is.na(x)]
  
  stats <- sapply(setdiff(sum_stats, "quantile"), function(fun){
    ifelse(length(x)==0, NA, do.call(fun, list(x, na.rm=na.rm)))
  })
  if(vname!="") names(stats) <- paste0(vname, "_", setdiff(sum_stats, "quantile"))
  
  if("quantile" %in% sum_stats){
    qq <- stats::quantile(x, probs = probs)
    names(qq) <- paste0("Q", probs)
    if(all(probs == c(0.25, 0.75))) names(qq) <- c("quart1", "quart3")
    if(vname!="") names(qq) <- paste0(vname, "_", names(qq))
    stats <- c(stats, qq)
  } 
  
  as.data.frame(t(stats))
}

get_hist <- function(x, vname="", nbins = 10){
  x <- x[!is.na(x)]
  
  binfreq <- nbinfreq <- c(rep(NA, nbins))
  if(length(x) > 1) {
    if(min(x) != max(x)){
      binfreq <- as.vector(graphics::hist(x,breaks=seq(from=min(x), to=max(x), 
                                             by=(max(x)-min(x))/nbins), plot = FALSE)$counts) 
      nbinfreq <- binfreq / length(x)
    }
  }
  
  names(binfreq) <- paste0("bin", 1:nbins)
  names(nbinfreq) <- paste0("nbin", 1:nbins)
  if(vname!="") names(binfreq) <- paste0(vname, "_", names(binfreq))
  if(vname!="") names(nbinfreq) <- paste0(vname, "_", names(nbinfreq))
  
  cbind(as.data.frame(t(binfreq)), as.data.frame(t(nbinfreq)))
}

################
################


attrWithOutliers <- function(ds,tgt, exclude = integer(0)) {
  
  numattrs <- as.numeric(which(sapply(ds,is.numeric)))
  exclude <- which(colnames(ds) %in% exclude)
  numattrs <- setdiff(numattrs,c(tgt, exclude))
  
  numAttrs <- 0
  if(length(numAttrs)>0)
    numAttrs <- sum(sapply(1:length(numattrs), function(i)
      length(grDevices::boxplot.stats(ds[,numattrs[i]])$out)))
  
  data.frame(attrWithOutliers = numAttrs)
}

########

typeOutliers <- function(ds,tgt) {
  
  stats <- grDevices::boxplot.stats(ds[,tgt])
  stats_dist <- stats$stats
  stats_out <- stats$out
  
  typeOutliersTgt <- "N"
  
  if(any(stats_out > stats_dist[5]) & any(stats_out < stats_dist[1])) {
    typeOutliersTgt <- "B"
  } else if(any(stats_out > stats_dist[5])) {
    typeOutliersTgt <- "U"
  } else if(any(stats_out < stats_dist[1])) {
    typeOutliersTgt <- "L"
  } else {}
  
  rm(stats)
  rm(stats_dist)
  rm(stats_out)
  
  data.frame(typeOutliersTgt = typeOutliersTgt)
  
}

#########

coefDet <- function(ds, formula, exclude = integer(0), nom=TRUE) {
  
  if(length(exclude)>0) ds <- ds[ , -exclude]
  
  if(nom) {
    m <- stats::lm(formula, data=ds)
  } else {
    numattrs <- as.numeric(which(sapply(ds,is.numeric)))
    new_ds <- ds[,numattrs]
    if(length(numattrs)>1) {
      m <- stats::lm(formula, data=new_ds)
    } else {
      m <- 0
    }
    
  }
  
  rsq <- 0
  if(is.numeric(m)) { } else {
    rsq <- summary(m)$r.squared
  }
  
  data.frame(coefDet = rsq)
  
}


#########

Cor_NumAttrs <- function(ds, tgt, exclude=integer(0),
                         sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), 
                         nbins=5, use="na.or.complete") {
  
  numattrs <- as.numeric(which(sapply(ds,is.numeric)))
  exclude <- which(colnames(ds) %in% exclude)
  numattrs <- setdiff(numattrs, c(tgt, exclude))
  
  v_cor <- stats::cor(ds[ , numattrs], use=use)
  v_cor[upper.tri(v_cor, diag = TRUE)] <- NA
  v_cor <- as.vector(v_cor)
  v_cor <- v_cor[!is.na(v_cor)]
  
  stats <- get_stats(v_cor, vname = "CorNumAttrs", sum_stats = sum_stats)
  if(nbins>0) hist <- get_hist(v_cor, nbins=nbins, vname = "CorNumAttrs")
  
  if(nbins>0) cbind(stats, hist)
  else stats
}


Cor_NumAttrs_Tgt <- function(ds, tgt, exclude=integer(0),
                             sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), 
                             nbins=5, use = "na.or.complete") {
  
  numattrs <- as.numeric(which(sapply(ds,is.numeric)))
  exclude <- which(colnames(ds) %in% exclude)
  numattrs <- setdiff(numattrs, c(tgt, exclude))
  
  tgt_cor <- stats::cor(ds[ ,tgt], ds[ ,numattrs], use=use)[1,] 
  tgt_cor <- tgt_cor[!is.na(tgt_cor)]

  stats <- get_stats(tgt_cor, vname = "CorTgt", sum_stats = sum_stats)
  if(nbins>0) hist <- get_hist(tgt_cor, nbins=nbins, vname = "CorTgt")
  
  if(nbins>0) cbind(stats, hist)
  else stats
  
}

#########

Moments_NumAttrs <- function(ds, tgt, exclude=integer(0), 
                             moments = c("geary", "kurtosis", "moment", "skewness"), 
                             m_labels = c(geary = "GKur", kurtosis = "PKur", 
                                          moment = "Mom", skewness = "Skew"),
                             sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), nbins=5, na.rm = TRUE) {
  
  numattrs <- as.numeric(which(sapply(ds,is.numeric)))
  numattrs <- setdiff(numattrs,c(tgt, exclude))
  
  stats <- list()
  for(m in moments){
    fun <- get(m, as.environment("package:moments"))
    moms <- apply(ds[ , numattrs], 2, fun, na.rm=na.rm)
    
    stats[[m]] <- get_stats(moms, vname = paste0(m_labels[m], "NumAttrs"), 
                       sum_stats = sum_stats)
    if(nbins>0) stats[[m]] <- cbind(stats[[m]], get_hist(moms, nbins=nbins, vname = paste0(m_labels[m], "NumAttrs")))
  }
  
  dplyr::bind_cols(stats)
}

Moments_Tgt <- function(ds, tgt,
                             moments = c("geary", "kurtosis", "moment", "skewness"), 
                             m_labels = c(geary = "GKur", kurtosis = "PKur", 
                                          moment = "Mom", skewness = "Skew")) {
  
  
  stats <- sapply(moments, function(m){
    fun <- get(m, as.environment("package:moments"))
    mom <- fun(ds[ , tgt])
  })
  names(stats) <- paste0(m_labels, "Tgt")
  
  as.data.frame(t(stats))
}

#########

Missing_NumAttrs <- function(ds, tgt, exclude=integer(0), 
                             nORp = 0.2,
                             sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), nbins=5) {
  
  numattrs <- as.numeric(which(sapply(ds,is.numeric)))
  numattrs <- setdiff(numattrs,c(tgt, exclude))
  
  percNAperCol <- apply(ds[, numattrs], 2,
                           function(y) length(which(is.na(y))) / length(y))
  percNAperRow <- apply(ds[, numattrs], 1,
                         function(y) length(which(is.na(y))) / length(y))
  
  stats <- c(get_stats(percNAperCol, vname = "percNAperCol"),
             get_stats(percNAperRow, vname = "percNAperRow"))
  if(nbins>0) stats <- c(stats, get_hist(percNAperCol, nbins=nbins, vname = "percNAperCol"),
                    get_hist(percNAperRow, nbins=nbins, vname = "percNAperRow"))
  
  dplyr::bind_cols(list(stats, 
    data.frame(percColWithTooManyNA = length(which(percNAperCol >= nORp)) / length(numattrs),
               percRowsWithTooManyNA = length(which(percNAperRow >= nORp)) / nrow(ds))))
}



#########


Dist_Locs <- function(stations, loc_id="station",
                      sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), nbins=5) {
  
  distm <- get_spatial_dist_mat(stations, loc_id)
  diag(distm) <- NA
  min_dists <- apply(distm, 1, function(x) min(x, na.rm=T))
  max_dists <- apply(distm, 1, function(x) max(x, na.rm=T))
  avg_dists <- apply(distm, 1, function(x) mean(x, na.rm=T))
  
  distm[upper.tri(distm, diag=T)] <- NA
  distm[upper.tri(distm, diag = TRUE)] <- NA
  distm <- as.vector(distm)
  distm <- distm[!is.na(distm)]
  
  dists <- list(MinDist = min_dists, MaxDist = max_dists, 
                AvgDist = avg_dists, Dist = distm)
  
  stats <- c()
  hist <- c()
  for(d in 1:length(dists)){
    stats <- c(stats, get_stats(dists[[d]], vname = names(dists)[d], sum_stats = sum_stats))
    hist <- get_hist(dists[[d]], nbins=nbins, vname = names(dists)[d])
  }
  sum_stats <- setdiff(sum_stats, c("max", "min"))
  nstats <- get_stats(distm / max(distm), 
                      vname = "Ndist", sum_stats = sum_stats)
 
  dplyr::bind_cols(list(stats, nstats, hist))
}


#########

ArimaStats <- function(ds, formula, time_id="time", loc_id="station",
                        sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), nbins=5,
                        ...){
  tgt <- as.character(as.list(formula)[[2]])

  arstats <- list()
  for(l in unique(ds[[loc_id]])){
    x <- ds[which(ds[[loc_id]] == l), c(time_id, loc_id, tgt)]
    x[[time_id]] <- as.POSIXct(x[[time_id]])
    
    ts <- xts::xts(x[[tgt]], order.by = x[[time_id]])
    period <- as.numeric(xts::periodicity(ts)$difftime)
    
    colnames(x) <- c("time", "station", "value")
    x <- x %>% tidyr::complete(time = seq.POSIXt(min(.data$time), max(.data$time), by=period))
    ts <- xts::xts(x[[tgt]], order.by = x[[time_id]])
      
    m <- forecast::auto.arima(ts)
    # m <- forecast::auto.arima(ts, ...)
    arma_labels <- c("ArimaAR", "ArimaMA", "ArimaSAR", "ArimaSMA", "ArimaT", "ArimaDiff", "ArimaSDiff")
    arstats[[l]] <- c(m$arma, m$aic, length(ts), length(which(is.na(ts))), length(which(is.na(ts))) / length(ts))
    names(arstats[[l]]) <- c(arma_labels, "aic", "TSLenPerLoc", "MissingTSPerLoc", "MissingTSPerLocPerc")  
    
  }
  arstats <- dplyr::bind_rows(lapply(arstats, function(x) as.data.frame(t(x))))
  
  stats <- dplyr::bind_cols(lapply(colnames(arstats), function(c) get_stats(arstats[ ,c], 
                                                                  sum_stats = sum_stats, vname=c)))
  if(nbins>0) hist  <- dplyr::bind_cols(lapply(colnames(arstats), function(c) get_hist(arstats[ ,c], 
                                                                 nbins = nbins, vname=c)))
  
  if(nbins>0) cbind(stats, hist)
  else stats
}

#########

MoranStats <- function(ds, stations, formula, time_id="time", loc_id="station",
                  p.thr = 0.05, sum_stats=c("min", "median", "quantile", "mean", "max", "sd"), nbins=5) {
  
  tgt <- which(colnames(ds) == as.list(formula)[[2]])
  
  distm <- get_spatial_dist_mat(stations, loc_id)
  
  morans <- list()
  for(t in unique(ds[[time_id]])){
    
    x <- ds[which(ds[[time_id]] == t), ]
    x <- x[order(x[[loc_id]]), ]
    ss <- x[ , tgt]
    if(nrow(x) > 1 & length(unique(ss))>1){
      dd <- distm[paste0("SITE_", x[ ,loc_id]),paste0("SITE_", x[ ,loc_id])]
      dd <- 1/dd
      diag(dd) <- 0
      dd[is.infinite(dd)] <- 0
      
      m <- list(observed=NA, p.value=NA)
      err <- try({ m <- ape::Moran.I(ss, dd, scaled = FALSE) })
      # if(class(err) == "try-error") print(paste("error at", t))
      # if(class(err) == "character") print(paste("warning at", t))
      
      morans[[as.character(t)]] <- data.frame(MoranObs = m$observed, 
                                              MoranObsSig = ifelse(m$p.value < p.thr, m$observed, NA),
                                              MoranPVal = ifelse(m$p.value < p.thr, 1, 0))
    }
  }
  morans <- dplyr::bind_rows(morans)
  
  stats <- dplyr::bind_cols(lapply(c("MoranObs", "MoranObsSig"), function(c) get_stats(morans[ ,c], 
                                                                            sum_stats = sum_stats, vname=c)))
  if(nbins>0) hist  <- dplyr::bind_cols(lapply(c("MoranObs", "MoranObsSig"), function(c) get_hist(morans[ ,c], 
                                                                           nbins = nbins, vname=c)))
  
  sig <- data.frame(MoranSigPval = sum(morans$MoranPVal, na.rm=T),
                    MoranSigPvalPerc = sum(morans$MoranPVal, na.rm=T) / nrow(morans),
                    MoranObsSigPos = sum(length(which(morans$MoranObsSig > 0))),
                    MoranObsSigPosPerc = sum(length(which(morans$MoranObsSig > 0))) / nrow(morans),
                    MoranObsSigNeg = sum(length(which(morans$MoranObsSig < 0))),
                    MoranObsSigNegPerc = sum(length(which(morans$MoranObsSig < 0)/ nrow(morans))))
  
  if(nbins>0) dplyr::bind_cols(list(stats, hist, sig)) 
  else dplyr::bind_cols(list(stats, sig))
}

#########

getMetaFeatures <- function(ds, formula, stations, time_id = "time", loc_id = "station",
                         nbins = 5){
  
  tgt <- which(colnames(ds) == as.list(formula)[[2]])
  coords <- which(colnames(ds) %in% c(time_id, loc_id))
  
  
  # counts
  
  counts <- data.frame(nCases = nrow(ds),
              nTimestamps = length(unique(ds$time)),
              nStations = length(unique(ds$station)),
              nAttr = ncol(ds)-3,
              numVars = length(as.numeric(which(sapply(ds[,c(-tgt, -coords)],is.numeric)))))
  
  # extreme stats
  
  ph <- uba::phi.control(ds[ ,tgt], method="extremes", coef=1.5)
  ls <- uba::loss.control(ds[ ,tgt])
  phi <- uba::phi(y = ds[ ,tgt], phi.parms = ph)
  
  extr_stats <- data.frame(nBoxplotExtremes = length(grDevices::boxplot.stats(ds[,tgt])$out),
    extremeRatio = length(grDevices::boxplot.stats(ds[,tgt])$out)/counts[["nCases"]],
    extremeRatio_0.9 = as.numeric( length(which(phi>=0.9)) / counts[["nCases"]] ),
    extremeRatio_0.5 = as.numeric( length(which(phi>=0.5)) / counts[["nCases"]] ),
    ratioCasesAttr = counts[["nCases"]] / counts[["nAttr"]],
    ratioTimesLocs = counts[["nTimestamps"]] / counts[["nStations"]],
    attrWithOutliers = attrWithOutliers(ds,tgt)[1,1] / counts[["nAttr"]])
  
  # spatio-temporal stats
  
  st_stats <- data.frame(Avail = 100*counts[["nCases"]] / (counts[["nTimestamps"]]*counts[["nStations"]]))
  
  time_stats <- ArimaStats(ds, formula, time_id, loc_id, nbins = nbins)
  space_stats <- cbind(Dist_Locs(stations, loc_id, nbins = nbins), 
                       MoranStats(ds, stations, formula, time_id, loc_id, nbins = nbins) )#,
                       #MoranStats2(ds, stations, formula, time_id, loc_id, nbins = nbins))
  
  st_stats <- dplyr::bind_cols(list(st_stats, time_stats, space_stats))
  
  # more complex stats
  
  if(length(coords) > 0){
    x <- ds[ , -coords] 
    x_tgt <- which(colnames(x) == as.list(formula)[[2]])
  }else{ 
    x <- ds
    x_tgt <- tgt
  }
  
  tgt_stats  <- bind_cols(list(data.frame(coefVar = stats::sd(ds[,tgt]) / mean(ds[,tgt]),
              coefDetTgt = coefDet(x, formula, nom=FALSE),
              stationaryTgt = ifelse(stats::sd(ds[,tgt]) > mean(ds[,tgt]), 1, 0),
              typeOutliers = typeOutliers(x, x_tgt)),
              Cor_NumAttrs_Tgt(x, x_tgt, nbins = nbins), Moments_Tgt(x, x_tgt)))
  
  stats <- dplyr::bind_cols(lapply(c("Missing_NumAttrs", "Cor_NumAttrs", "Moments_NumAttrs"),
         function(fun) do.call(fun, list(x, x_tgt, nbins = nbins))))
  
  bind_cols(list(counts, st_stats, extr_stats, tgt_stats, stats))
}

