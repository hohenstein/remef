#' Find out whether a model includes an intercept
#' 
#' The function tests whether an intercept is present in a model.
#' 
#' @param model a model object.
#' 
#' @return logical. \code{TRUE} if the intercept is present,
#'         \code{FALSE} otherwise.
#' 
#' @export
#' @keywords internal

has_intercept <- 
function(model) 
{
  as.logical(attr(terms(model), "intercept"))    
}
