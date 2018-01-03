#' @title Residual sum of squares (RSS) for all time series in a group.
#' 
#' @description The trend function for each individual time series is estimated non-parametrically 
#' by the local linear estimate (as discussed in Fan and Gijbels (1996)). Then, the common trend 
#' function for the group is estimated by taking the average over the group. Next, the shift functions
#' are estimated as the differences from the individual trend functions and finally, the residual sum
#' of squares are calculated  using the original values, the common trend functions and the shifts.
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. Every 
#' data frame should have two columns with names 'distance' and 'oxygen'.
#'  
#' @param S A vector of integers showing which individuals are considered in the group.
#' 
#' @param bandwidth Denotes the order of the bandwidth that should be used in the estimation process. 
#' bandwidth = k will mean that the bandwidth is n^k.
#' 
#' @return A vector of length equal to the group-size, so that each element is the RSS for the
#' corresponding individual in the group.

calculateRSS <- function(paths,
                         S,
                         bandwidth){
  
  RSS = numeric(length(S))
  for (i in 1:length(S)){
    y = cbind(paths[[S[i]]]$oxygen,log(paths[[S[i]]]$distance))
    n = dim(y)[1]
    ndx = (1:n)/n
    trend = array(0,dim=c(length(S),n,2))
    for (j in 1:length(S)){
      z = paths[[S[j]]]$oxygen
      trend[j,,1] = EstTrend(z,ndx,bandwidth)
      z = log(paths[[S[j]]]$distance)
      trend[j,,2] = EstTrend(z,ndx,bandwidth)
    }
    avgtrend = apply(trend,c(2,3),mean)
    ck = colMeans(trend[i,,]-avgtrend)
    RSS[i] = 0
    for (k in 1:n){
      RSS[i] = RSS[i] + sum((y[k,]-avgtrend[k,]-ck)^2)
    }
  }

  return(RSS)

}
