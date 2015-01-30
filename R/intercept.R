#' Get the intercept of a model
#' 
#' The function returns the value of the beta coeffcient of the model's 
#' intercept.
#' 
#' @param model a model object.
#' 
#' @return a length-one numeric vector.
#' 
#' @details \code{model} is an object returned by a model fitting 
#'   function (e.g., \code{\link{lm}}, \code{\link{aov}}).
#'   
#'   If the model was fitted without an intercept, the function
#'   returns \code{0} and generates a warning message.
#'   
#' @keywords internal
#' @examples
#' require(utils)
#'  
#' data(iris)
#' fit <- lm(Sepal.Length ~ Sepal.Width * Petal.Length * Petal.Width, data = iris)
#' summary(fit)
#' remef:::intercept(fit)

intercept <- 
function(model)
{
  if ( !has_intercept(model) ) {
    warning("The model was fitted without an intercept.")
    int <- 0
  } else {
    moma <- model.matrix(model)
    asgn <- attr(moma, "assign")
    idx <- asgn == 0L
    if ( is(model, "merMod") )
      coefs <- fixef(model)
    else
      coefs <- coef(model)
    int <- unname(coefs[idx])
  }
  int
}
