#' Calculate random effects based on a subset of variance terms
#'
#' This function calculates random effects of a mixed model based
#' on a subset of the variance terms.
#' 
#' @param ran a named list of integer vectors representing
#'   the random effects that should be calculated.
#' @param model a model object.
#' 
#' @return a numeric vector.
#' 
#' @details The names of the list \code{ran} correspond to the names
#'   of the list returned by \code{ranef_labels(model)}.
#'   
#'   The integer values within each vector correpsond to the position
#'   of the name of the random effect.
#'   
#'   If \code{ran} is \code{NULL}, the function will return a vector
#'   of zeroes.
#'                        
#' @seealso \code{\link{ranef}} for extracting the modes of the
#'   random effects.
#'   
#' @export
#' @keywords internal

calc_ranef <-
function(model, ran)
{
  if ( length(ran) == 0L )
    return( numeric(length = getME(model, "n")) )
  if ( is.list(ran) && length(ran) > 0L ) {
    stopifnot(is.numeric(unlist(ran)))
    ran <- lapply(ran, unique)
    ran_valid_idx <- vapply(ran, length, integer(1)) > 0L
    ran <- ran[ran_valid_idx]
  }
  ran_labs <- ranef_labels(model)
  rf_labs <- names(ran_labs)
  if ( is.null(names(ran)) ||
         any( !names(ran) %in% rf_labs ) ||
         any(duplicated(names(ran)))
  )
    stop("The list 'ran' has invalid names.")
  # number of random effects for each random factor
  num_re <- vapply(ran_labs, length, FUN.VALUE = integer(1L),
                   USE.NAMES = FALSE)
  num_re_previous <- c(0L, cumsum(head(num_re, -1L)))
  rf_idx <- match(names(ran), names(ran_labs))
  ran_model <- ranef(model)
  # subset of random effects corresponding to 'ran'
  ran_model_sub <- lapply(mapply("[", ran_model[rf_idx], ran,
                                 SIMPLIFY = FALSE),
                          unlist)
  # transform the list of 'Zt' matrices into a list of lists
  re_list_idx <- mapply(function(x, y) seq_len(x) + y,
                        num_re, num_re_previous,
                        SIMPLIFY = FALSE)
  re_list_flat <- getME(model, "Ztlist")
  re_list <- lapply(re_list_idx, function(x) re_list_flat[x])
  # extract subsets
  re_list_sub <- mapply("[", re_list[rf_idx], ran,
                        SIMPLIFY = FALSE)
  ran_model_sub <- mapply("[", ran_model[rf_idx], ran,
                          SIMPLIFY = FALSE)
  # here the calculation takes place
  ranef_prods <- mapply(function(mats, vecs)
    mapply(function(v, m) as.vector(v %*% m),
           vecs, mats),
    re_list_sub, ran_model_sub,
    SIMPLIFY = FALSE)
  ranef_sums <- rowSums(do.call(cbind, ranef_prods))
  ranef_sums
}
