#' @title Splitting-Coalescence (SC) algorithm.
#' 
#' @description Performs the iterative clustering algorithm on the archaeological time series data.
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. There 
#' should be two columns with names 'distance' and 'oxygen'.
#' 
#' @param bandwidth Denotes the order of the bandwidth that should be used in the splitting-coalescence 
#' (SC) clustering algorithm. A value k will mean that the bandwidth used in the algorithm is n^k.
#' 
#' @return A list of vectors where each vector gives the indexes of the individuals to be assigned
#' in the same cluster.

SCalgo <- function(paths,
                   bandwidth){
  
  p = length(paths)
  U = 1:p
  S = list()
  i = 1
  while (length(U)>0){
    out = iteration(paths,U,bandwidth)
    S[[i]] = out$S1
    i = i+1
    U = out$U
  }
  
  return(S)
  
}