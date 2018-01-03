#' @title Splitting-Coalescence-Estimation Method (SCEM) for archaeological time series.
#' 
#' @description Performs the clustering algorithm SCEM on the bivariate time series data - where
#' one series is for the distance from the cementum-enamel junction, and the other series is for
#' the value of the oxygen-18 isotope at that distance - and returns the class-assignments and 
#' birth seasobality estimates for all the individuals. 
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. There 
#' should be two columns with names 'distance' and 'oxygen'.
#' 
#' @param bandwidth Denotes the order of the bandwidth that should be used in the splitting-coalescence 
#' (SC) clustering algorithm. A value k will mean that the bandwidth used in the algorithm is n^k.
#' 
#' @return A list that contains a data frame and a list of vectors. The data frame has the individual 
#' information (ID, species, number of observations in the time series), cluster assignment, estimated 
#' period, delay and the birth seasonality estimate for every individual. The list of vectors gives the
#' groups formed by the clustering algorithm in the method.

SCEM <- function(paths,
                 bandwidth){
  
  cluster = SCalgo(paths,bandwidth = bandwidth)
  groups = cluster
  cosine = makeFits(paths)
  period = cosine$X
  gnum = length(groups)
  birthseason = numeric(gnum)
  for (k in 1:gnum){
    S = groups[[k]]
    for (jj in 1:length(S)){
      paths[[S[jj]]]$zval = paths[[S[jj]]]$distance/period[S[jj]]
    }
    fulldata = do.call(rbind,paths[S])
    fulldata = fulldata[order(fulldata$zval),]
    ddd = aggregate(oxygen ~ zval,fulldata,FUN = mean)
    xx = (ddd$zval)
    yy = ddd$oxygen
    zz = xx[sort.int(yy,decreasing = T,index.return = T)$ix[1:5]]
    if (max(dist(zz))>0.5){
      zz[zz<0.5] = 1+zz[zz<0.5]
    }
    birthseason[k] = mean(zz)%%1
  }
  
  index <- numeric(length(paths))
  season <- numeric(length(paths))
  species <- unique(do.call(rbind,paths)[,"ID"])
  results <- data.frame(ID = species)
  results$Length = as.numeric(lapply(paths,nrow))
  dd = matrix(nrow = length(paths),ncol = 2)
  for (i in 1:length(groups)){
    dd[groups[[i]],1] = as.numeric(i)
    dd[groups[[i]],2] = as.numeric(round(birthseason[i],3))
  }
  
  results$Cluster = dd[,1]
  results$x0 = dd[,2]*period
  results$X = period
  results$Season = dd[,2]
  
  return(list(results = results,groups = groups))
}
