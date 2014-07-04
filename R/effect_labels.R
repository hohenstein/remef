#' Get the labels of the fixed effects of the model
#' 
#' The function returns the labels of coefficients of the model.
#' Note: These labels are not identical to the term labels.
#' 
#' @param model a model object.
#' 
#' @return a character vector. 
#' 
#' @export
#' @keywords internal

effect_labels <- 
function(model) 
{
  colnames(model.matrix(model))    
}