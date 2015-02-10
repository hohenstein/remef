#' Find associated effects
#'
#' This function finds the labels of coefficients associated with
#' interactions or main effects in the results of various model
#' fitting functions.
#' 
#' @param model a model object.
#' @param ef a character string (a coefficent label).
#' @param order a character string vector giving the order of the
#'   terms to be returned.
#'   One of \code{"higher"} (default), \code{"lower"},
#'   or \code{"all"}, can be abbreviated.
#' @param include.base logical. Should the labels in \code{ef}
#'   be included in the output?
#'
#' @return \code{asef} returns vectors of character strings
#'   (coefficent labels).   
#'                        
#' @details \code{model} is an object returned by a model fitting
#'   function (e.g., \code{\link{lm}}, \code{\link{aov}}).
#'   
#'   The character string passed to \code{ef} must not include the 
#'   intercept (\code{"(Intercept)"}) since it cannot be part of an
#'   interaction. All non-intercept labels of any order are allowed.
#'   
#'   If \code{order} is \code{"higher"}, all terms of a higher order
#'   are returned. If \code{order} is \code{"lower"}, all terms of a 
#'   lower order are returned. If \code{order} is \code{"all"} then
#'   both terms of lower and higher order are returned. For example,
#'   if a model fit has the terms \code{A}, \code{B}, \code{C},
#'   \code{A:B}, \code{B:C}, \code{A:C}, and \code{A:B:C}, terms of 
#'   higher order relative to \code{A}  are \code{A:B}, \code{A:C}, 
#'   and \code{A:B:C}. Terms of lower order relative to \code{A:B} are
#'   \code{A} and \code{B}.
#'   
#'   If \code{include.base} is \code{FALSE} (the default), the term
#'   label in \code{ef} is not part of the returned string vector.
#'   
#' @seealso \code{\link{terms}} for extracting model terms.
#' @export
#' @examples
#' require(utils)
#'  
#' data(iris)
#' fit <- lm(Sepal.Length ~ Sepal.Width * Petal.Length * Petal.Width, data = iris)
#' asef(fit, "Petal.Length:Petal.Width", order = "higher")
#' asef(fit, "Sepal.Width:Petal.Length:Petal.Width", order = "lower")
#' asef(fit, "Petal.Length:Petal.Width", order = "all")
#'  
#' data(warpbreaks)
#' fit2 <- lm(breaks ~ wool * tension, data = warpbreaks)
#' # Since wool is a factor, the coefficient has another name, "woolB"
#' asef(fit2, "woolB", order = "higher")
#' asef(fit2, "woolB", order = "higher", include.base = TRUE)

asef <- 
function(model, ef, order = c("higher", "lower", "all"), 
         include.base = FALSE)
{
  if ( length(ef) != 1L )
    stop("ef must include exactly one value.")
  order <- match.arg(order)
  moma <- model.matrix(model)
  mterms <- terms(model)
  full_labs <- effect_labels(model)
  asgn <- attr(moma, "assign")
  idx <- match(ef, full_labs)
  aidx <- asgn[idx]
  if ( !aidx )
    stop("Hierarchy of effects is not defined for the intercept.")
  var_mat <- attr(mterms, "factors")
  ef_order <- attr(mterms, "order")
  rel_matrix <- var_mat[ var_mat[ , aidx] > 0L, , drop = FALSE ]
  rel_cum <- colSums(rel_matrix)
  ef_order_base <- ef_order[aidx]
  lab_idx <- 
    which(switch (order,
                  higher = {
                    rel_cum == ef_order_base
                  },
                  lower = {
                    rel_cum == ef_order
                  },
                  all = {
                    rel_cum == ef_order_base | rel_cum == ef_order
                  }
    ))
  lab_list <- sapply(lab_idx, function(i) {
    if (ef_order[i] < ef_order_base) {
      ef_splitted <- strsplit(ef, ":", fixed = TRUE)[[1L]]
      combinations <- apply(combn(ef_splitted, ef_order[i]),
                            2L, paste, collapse = ":")
      full_labs[asgn == i][full_labs[asgn == i] %in% combinations]
    } else if (i == aidx) {
      if (include.base) ef
    } else {
      cand_labs <- full_labs[asgn == i]
      cand_splitted <- strsplit(cand_labs, ":", fixed = TRUE)
      ef_splitted <- strsplit(ef, ":", fixed = TRUE)[[1L]]
      matches <- do.call(cbind, 
                         lapply(cand_splitted, "%in%", x = ef_splitted))
      cand_labs[colSums(matches) == ef_order_base]
    }
  })
  unlist(lab_list, use.names = FALSE)
}
