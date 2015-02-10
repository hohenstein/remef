#' Extract the coefficient names of terms of a model.
#' 
#' The function returns the labels of the coefficients of the model
#' Note. The terms are not identical to the coefficients.
#' 
#' @param model a model object.
#' @param term a character vector of term labels.
#' @param as.list logical. Should the result be returned as a list?
#' 
#' @return either a character vector (if \code{as.list = FALSE})
#'   or a list with character vectors (otherwise) with one 
#'   list element for each element in \code{term}. 
#'
#' @details \code{model} is an object returned by a model fitting 
#'   function (e.g., \code{\link{lm}}, \code{\link{aov}}).
#'   
#'   The character string passed to \code{term} must not include the 
#'   intercept (\code{"(Intercept)"}).
#'   All non-intercept labels in any order are allowed.
#'   
#' @seealso \code{\link{coef2term}} for the inverse function.
#' @export
#' @examples
#' require(utils)
#'  
#' data(iris)
#' fit <- lm(Sepal.Length ~ poly(Sepal.Width, 3), data = iris)
#' summary(fit)
#' term2coef(fit, "poly(Sepal.Width, 3)")
#'  
#' data(warpbreaks)
#' fit2 <- lm(breaks ~ wool * tension, data = warpbreaks)
#' summary(fit2)
#' term2coef(fit2, c("wool", "tension"))
#' term2coef(fit2, c("wool", "tension"), as.list = TRUE)
#' term2coef(fit2, c("wool:tension"))

term2coef <- 
function(model, term, as.list = FALSE) 
{
  full_labs <- effect_labels(model)
  term_labs <- term_labels(model)
  moma <- model.matrix(model)
  asgn <- attr(moma, "assign")
  idx <- match(term, term_labs)
  new_labs <- lapply(idx, function(x) full_labs[asgn == x])
  if (!as.list)
    unlist(new_labs, use.names = FALSE)
  else
    setNames(new_labs, term)
}
