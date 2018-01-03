#' @title Estimates the trend function for a time series.
#' 
#' @description The trend function for each individual time series is estimated non-parametrically 
#' by the local linear estimate (as discussed in Fan and Gijbels (1996)). Detailed description can 
#' be found in the supplementary file of the paper.
#' 
#' @param y A vector of time series observations.
#' 
#' @param time A vector of time points where the value of the trend needs to be estimated.
#' 
#' @param bandwidth Denotes the order of the bandwidth that should be used in the estimation process. 
#' bandwidth = k will mean that the bandwidth is n^k.
#' 
#' @return A vector of estimated values for the trend function at the given time-points.

EstTrend <- function(y,
                     time,
                     bandwidth){
  
  n = length(y)
  bn = n^{bandwidth}
  x = 1:n
  out = numeric(length(time))
  for (i in 1:length(time)){
    t = time[i]
    S0 = sum(kernel((x/n-t)/bn))
    S1 = sum((t-x/n)*kernel((x/n-t)/bn))
    S2 = sum((t-x/n)^2*kernel((x/n-t)/bn))
    w = kernel((x/n-t)/bn)*(S2-(t-x/n)*S1)/(S2*S0-S1^2)
    out[i] = sum(y*w)
  }
  
  return(out)

}
