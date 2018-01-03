#' @title Prepare results for cosine model fit.
#' 
#' @description Performs the nonlinear least squares (NLS) regression method for the cosine 
#' model, with the proposed initialization for all the parameters. It fits the NLS method 
#' as required, and then computes different quantities for the birth seasonality estimates
#' corresponding to different individuals.
#' 
#' @param paths A list of data frames, where each frame contains the data for one individual. Every 
#' data frame should have two columns with names 'distance' and 'oxygen'.
#' 
#' @return A data frame with the estimated parameters in the model, birth seasonality estimate, 
#' predicted/observed minimum/maximum for the oxygen isotope variable, mean squared error
#' and Pearson's R^2 corresponding to the model fit for every individual.

makeFits = function(paths) {
  
  fits = c()
  for(i in 1:length(paths)) {
    data = paths[[i]]
    curve = sineFit(data)
    fit = convertParameters(curve)
    fit$predictedMin = fit$intercept - abs(fit$amplitude)
    fit$predictedMax = fit$intercept + abs(fit$amplitude)
    fit$observedMin = min(data$oxygen)
    fit$observedMax = max(data$oxygen)
    fit$MSE = mean((predict(curve) - data$oxygen)^2)
    fit$PearsonCorrelation = cor(predict(curve),paths[[i]]$oxygen,method = "pearson")
    fits = rbind(fits, as.numeric(fit))
  }
  fits = data.frame(fits)
  colnames(fits) = c("amplitude","intercept","x0","X","birth","predictedMin",
                     "predictedMax","observedMin","observedMax","MSE","Pearson")
  
  return(fits)
  
}
