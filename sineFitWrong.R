#' @title Cosine model fitting with given initialization for two parameters.
#' 
#' @description Performs the updated nonlinear least squares (NLS) regression method for the cosine 
#' model proposed by Balasse et al. The method calculates with the proposed initial values at first, 
#' and then fits the NLS method as required.
#' 
#' @param data A data frame that contains the data for one individual. There should be two columns 
#' with names 'distance' and 'oxygen'.
#' 
#' @param amplitude A number, corresponding to a given value for the amplitude parameter.
#' 
#' @param intercept A number, corresponding to a given value for the intercept parameter.
#' 
#' @return A fitted model object from the nls function in R.

sineFitWrong = function(data,
                        amplitude,
                        intercept) {
  
  frequency = 2*pi/max(data$distance)
  model = lm(oxygen ~ sin(frequency*distance) + cos(frequency*distance),data = data)
  a = model$coef[2]
  b = model$coef[3]
  phase = -atan(a/b)
  
  start = list(intercept = intercept, amplitude = amplitude,phase = phase,frequency = frequency)
  curve = nls(oxygen ~ intercept + amplitude*cos(frequency*distance+phase),
              data = data,start = start,control = nls.control(warnOnly=TRUE))
  
  return(curve)
  
}
