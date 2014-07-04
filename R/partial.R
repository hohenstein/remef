#' Calculate partial effects
#'
#' \code{partial} is used to calculate partial effects of a (generalized) 
#'  linear mixed-effects model fitted by a function of the \code{lme4}
#'  package.
#' 
#' @param model a model object.
#' @param fix a character of coefficient labels (fixed effects).
#' @param ran a named list of integer vectors representing
#'   random effects.
#' @param grouping logical. Should terms associated with the ones
#'   in \code{fix} also be selected? The functionality of this
#'   parameter depends on the value of \code{keep}. 
#'   The details are given under 'Details'.
#' @param keep.intercept logical. Should the intercept be kept
#'   (default) or be removed?
#' @param inverse logical. If \code{TRUE}, the inverse of the link
#'   function is applied (after the effects are removed).
#' @param keep logical. Should the specified effects be removed
#'   (default) or be kept?
#'   
#' @return \code{partial} returns a numeric vector. The length equals
#'   the number of observations in the model.
#'                        
#' @details \code{model} is an object returned by a model fitting 
#'   function of the \pkg{lme4} package.
#'   The function works with both linear mixed-effects models
#'   (returned by the function \code{lmer}) and generalized linear
#'   mixed-effects models (\code{glmer}).
#'     
#'   The model formula for a (generalised) linear mixed-effects model
#'   can be written as
#'   
#'   \deqn{y = X \beta + Z b + \epsilon}{
#'         y = X * beta + Z * b + epsilon}
#'   
#'   where \eqn{y} is the vector including the response values,
#'   \eqn{X} is the (fixed-effects) design matrix, \eqn{\beta}{beta} 
#'   is the vector of population coefficients (fixed effecs), 
#'   \eqn{Z} is the random-effects design matrix, \eqn{b} is
#'   a vector of random effects, and \eqn{\epsilon}{epsilon} is a 
#'   vector of random error terms (residuals).
#'   
#'   The response vector can be represented as a sum of fixed effects,
#'   random effects, and residuals. In order to construct partial effects,
#'   the function \code{partial} either "removes" effects (if 
#'   \code{keep = FALSE}, see also \code{\link{remef}}) or "keeps"
#'   effects (if \code{keep = TRUE}, see also \code{\link{keepef}}).
#'   In both cases, the function does keep the residuals. Hence,
#'   the function does not return fitted values.
#'   
#'   The function removes a subset of fixed and random effects from
#'   the response variable. If all effects are removed, only the residuals
#'   remain.
#'   
#'   The valid strings for the parameter \code{fix} correspond to the
#'   names of the coeffcients in the model summary output. Note that
#'   the name of a coefficient is not necessarily identical to the name
#'   of a predictor variable (e.g., factor variables). Furthermore,
#'   a single variable in a model formula can result in multiple 
#'   model coefficients (e.g., polynomial contrasts).
#'   The character string passed to \code{fix} must not include the 
#'   intercept (\code{"(Intercept)"}). Use the parameter 
#'   \code{keep.intercept} instead.
#'   
#'   The parameter \code{ran} is used to specify random effects.
#'   The names of the list passed to \code{ran} correspond to the names
#'   of the random factors in the model summary.
#'   Alternatively, the string \code{"all"} can pe passed to \code{ran}
#'   which results in the selection of all random effects.
#'   
#'   If \code{grouping} is \code{FALSE}, only the specified fixed effects in
#'   \code{fix} are used. 
#'   If \code{grouping} is \code{TRUE}, associated effects of lower 
#'   \emph{or} higher order are used too. The actual behaviour of the 
#'   function depends on the \code{keep} parameter: If \code{keep = FALSE} 
#'   (effects are "removed"), \code{grouping} specifies whether effects of 
#'   higher order should be selected and removed as well. Otherwise, if 
#'   \code{group = TRUE} (effects are "kept"), \code{grouping} specifies 
#'   whether effects of higher order should be selected and removed as well. 
#'   For example,
#'   if a model fit has the terms \code{A}, \code{B}, \code{C},
#'   \code{A:B}, \code{B:C}, \code{A:C}, and \code{A:B:C}, terms of 
#'   higher order relative to \code{A}  are \code{A:B}, \code{A:C}, 
#'   and \code{A:B:C}. Terms of lower order relative to \code{A:B} are
#'   \code{A} and \code{B}.
#'   
#'   The parameter \code{inverse} is only important for \emph{generalised}
#'   linear mixed-effects models (in contrast to linear mixed-effects 
#'   models). In the former type of statistical models, a link function
#'   is specified to transform the response variable.
#'   If \code{partial} is used with such a type of model, the resulting
#'   values will be in the metric of the variable \emph{after}. For example,
#'   if you have a
#'   binary response variable and fit a model with the binomial distribution
#'   using the logit link funtion, the use of \code{partial} will return
#'   logits unless \code{inverse = TRUE} in which case logits will be 
#'   transformed to proportions after the removal of effects. 
#'            
#' @seealso \code{\link{remef}} for a wrapper with \code{keep = FALSE} 
#'   and \code{\link{keepef}} for a wrapper with \code{keep = TRUE}.
#' @export
#' @examples
#' # TODO: Some examples

partial <- 
function(model, fix = NULL, ran = NULL, grouping = FALSE,
         keep.intercept = TRUE, inverse = FALSE, keep = FALSE) 
{
  if ( !is(model, "merMod") ) 
    stop("This class is not supported yet.")
  
  mclass <- class(model)
  if ( keep || mclass == "glmerMod" ) { 	
    DV <- residuals(model)
  } else
    DV <- getME(model, "y")  
  
  full_labs <- effect_labels(model)
  if ( has_intercept(model) )
    full_labs <- full_labs[-1L]
  
  # Part 1) fixed effects
  eidx <- is.na(match(fix, full_labs))
  if ( any(eidx) )
    stop("The following effects are not present in the model:
          \t", paste(fix[eidx], collapse = ", "))
    
  if ( grouping && !is.null(fix) ) {
    order <- if (keep) "lower" else "higher"
    fix <- unique(unlist(sapply(fix, asef, model = model, order = order,
                                include.base = TRUE, simplify = FALSE)))    
  }
  
  if ( !keep && mclass == "glmerMod" )
    # complement (the 'full_labs' without 'fix')
    fix <- setdiff(full_labs, fix)
  
  moma <- model.matrix(model)
  fixef_prods <- as.vector(moma[ , fix, drop = FALSE] %*% 
                             fixef(model)[fix])
  
  
  if( keep || mclass == "glmerMod" ) {
    int <- if ( keep.intercept ) intercept(model) else 0    
    DV <- DV + fixef_prods + int
  } else {
    int <- if ( !keep.intercept ) intercept(model) else 0    
    DV <- DV - fixef_prods - int
  }
  
  # Part 2) random effects
  ran_labs <- ranef_labels(model)
  
  if ( identical(ran, "all") )
    ran <- lapply(ran_labs, seq_along)
  
  if ( !keep && mclass == "glmerMod" && length(ran) > 0L ) {
    # complement (the 'ran_labs' without 'ran')
    rf_labs <- names(ran_labs)
    
    ran_inv <- lapply(rf_labs, function(x) {
      idx <- match(x, names(ran))
      re_seq <- seq_along(ran_labs[[x]])
      if (is.na(idx)) {
        re_seq
      } else {
        setdiff(re_seq, ran[[x]])
      }})
    names(ran_inv) <- rf_labs
    ran <- ran_inv      
  }
  
  ranef_sums <- calc_ranef(ran, model)  
  
  if( keep || mclass == "glmerMod" ) {
    DV <- DV + ranef_sums
  } else {
    DV <- DV - ranef_sums
  } 
  
  # apply inverse function
  if ( inverse ) {
    fam <- family(model)
    DV <- fam$linkinv(DV)
  }
  unname(DV)
}