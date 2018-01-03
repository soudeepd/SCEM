#' @title Bayesian Information Criterion (BIC) for a partition.
#' 
#' @description This function calculates an extended version of BIC, which is computed using
#' a particular weighted average of the total residual sum of squares and the number of clusters.
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. Every 
#' data frame should have two columns with names 'distance' and 'oxygen'.
#'  
#' @param partition A list of vectors. Each element in the list is a vector of integers, corresponding
#' to individuals considered in one group.
#' 
#' @param bandwidth Denotes the order of the bandwidth that should be used in the estimation process. 
#' bandwidth = k will mean that the bandwidth is n^k.
#' 
#' @return Value of the extended BIC function for the partition.

EBIC <- function(paths,
                 partition,
                 bandwidth){
  
  q = length(partition)
  p = length(paths)
  size = as.numeric(lapply(paths,nrow))
  n = mean(size)
  bn = n^{bandwidth}
  total = 0
  for (i in 1:q){
    RSS = sum(calculateRSS(paths,partition[[i]],bandwidth))
    total = total+RSS
  }
  out = n*p*log(total/(n*p))+q*log(n*bn)*(1/bn-1)
  
  return(out)

}
