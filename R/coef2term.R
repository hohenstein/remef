#' Extract the term names of coefficients of a model.
#' 
#' The function returns the labels of the terms of the model
#' Note. The terms are not identical to the coefficients.
#' 
#' @param model a model object.
#' @param ef a character vector of coefficient labels.
#' 
#' @return a character vector of term labels (one element for each
#'    element in \code{ef}). 
#'
#' @details \code{model} is an object returned by a model fitting 
#'   function (e.g., \code{\link{lm}}, \code{\link{aov}}).
#'   
#'   The character string passed to \code{ef} must not include the 
#'   intercept (\code{"(Intercept)"}).
#'   All non-intercept labels of any order are allowed.
#'   
#' @seealso \code{\link{term2coef}} for the inverse function.
#' @export
#' @examples
#' require(utils)
#'  
#' data(iris)
#' fit <- lm(Sepal.Length ~ poly(Sepal.Width, 3), data = iris)
#' summary(fit)
#' coef2term(fit, "poly(Sepal.Width, 3)2")
#'  
#' data(warpbreaks)
#' fit2 <- lm(breaks ~ wool * tension, data = warpbreaks)
#' summary(fit2)
#' coef2term(fit2, c("woolB", "tensionH"))
#' coef2term(fit2, c("woolB:tensionM"))

coef2term <-
function(model, ef)
{
  full_labs <- effect_labels(model)
  term_labs <- term_labels(model)
  moma <- model.matrix(model)
  asgn <- attr(moma, "assign")
  aidx <- match(ef, full_labs)
  term_labs[asgn[aidx]]
}
