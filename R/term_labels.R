#' Get the labels of the terms of a model.
#' 
#' The function returns the labels of the terms of the model
#' Note. The terms are not identical top the coefficients.
#' 
#' @param model a model object.
#' 
#' @return a character vector. 
#' 
#' @export
#' @keywords internal

term_labels <- 
function(model) 
{
  labels(terms(model))  
}