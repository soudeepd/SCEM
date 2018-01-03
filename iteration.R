#' @title Iteration step for the Splitting-Coalescence-Estimation Method (SCEM).
#' 
#' @description This function performs the iterative step in the proposed method. 
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. Every 
#' data frame should have two columns with names 'distance' and 'oxygen'.
#'  
#' @param U A list of vectors. Each element in the list is a vector of integers, corresponding
#' to individuals considered in one group.
#' 
#' @param bandwidth Denotes the order of the bandwidth that should be used in the estimation process. 
#' bandwidth = k will mean that the bandwidth is n^k.
#' 
#' @return A list of two sets - the first one is a cluster and the second one is the set of
#' individuals to be used in the next iteration.

iteration <- function(paths,
                      U,
                      bandwidth){
  
  V = numeric()
  backup = U
  while (length(U)>0){
    RSS = calculateRSS(paths,U,bandwidth)
    M = which.max(RSS)
    V = c(V,U[M])
    U = setdiff(U,U[M])
  }
  
  ebicvalues = numeric(length(backup))
  ebicvalues[1] = EBIC(paths,list(backup),bandwidth)
  if (length(backup)>1){
    for (i in 2:length(backup)){
      partition = list()
      for (j in 1:(i-1)){
        partition[[j]] = V[j]
      }
      partition[[i]] = V[i:length(backup)]
      ebicvalues[i] = EBIC(paths,partition,bandwidth)
    }
  }  
  
  L = which.min(ebicvalues)
  S1 = V[L:length(backup)]
  U2 = setdiff(V,S1)
  out = list(S1,U2)
  names(out) = c("S1","U")
  
  return(out)
  
}
