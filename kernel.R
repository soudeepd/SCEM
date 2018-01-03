#' @title Epanechnikov kernel.
#' 
#' @description Calculates the value of the Epanechnikov kernel function for any vector.
#' 
#' @param v A vector of real numbers.
#' 
#' @return A vector of the calculated kernel values for the input vector.

kernel <- function(v){
  
  return(0.75*(1-v^2)*(abs(v)<1))

}