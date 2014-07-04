#' Get the labels of the random effects of the model.
#' 
#' The function returns the labels of coefficients of the model.
#' Note. These labels are not identical to the term labels.
#' 
#' @param model a model object.
#' 
#' @return a named list of character vectors. 
#' 
#' @export
#' @keywords internal

ranef_labels <- 
function(model) 
{
  getME(model, "cnms")
}