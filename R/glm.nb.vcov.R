#' Extended variance matrix
#'
#' An extension to the \code{\link[stats]{vcov}} function mainly to
#' cover the additional parameter involved in negative binomial models.
#' (Currently the same as \code{\link[stats]{vcov}} apart from negative binomial models.)
#'
#' @param object A fitted mdel objeds
#' @param ... currently ignored
#'
#' @return An extended variance matrix including parameters addition to the regression coefficients
#' @export
#'
#' @examples
#' fm <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), quine)
#' Sigma <- vcovx(fm)
vcovx <- function(object, ...) {
  UseMethod("vcovx")
}

#' @rdname vcovx
#' @export
vcovx.default <- function(object, ...) {
  stats::vcov(object, ...)
}

#' @rdname vcovx
#' @export
vcovx.negbin <- function(object, ...) {

  ## given a model fitted by glm.nb in MASS, this function returns a
  ## variance covariance matrix for the regression coefficients and
  ## dispersion parameter, without assuming independence between these
  ## note that the model must have been fitted with x=TRUE argument so
  ## that design matrix is available

  ## formulae based on p23-p24 of
  ## http://pointer.esalq.usp.br/departamentos/lce/arquivos/aulas/2011/LCE5868/OverdispersionBook.pdf
  ## and
  ## http://www.math.mcgill.ca/~dstephens/523/Papers/Lawless-1987-CJS.pdf

  ## Author: Jonathan Bartlett,
  ## https://stats.stackexchange.com/questions/221648/negative-binomial-regression-in-r-allowing-for-correlation-between-dispersion

  if(!("x" %in% names(object))) {
    Call <- substitute(update(MOD, x = TRUE), list(MOD = substitute(object)))
    object <- eval(Call, parent.frame())
  }

  k <- object$theta
  ## p is number of regression coefficients
  p <- nrow(vcov(object))

  ## construct observed information matrix
  obsInfo <- matrix(0, p + 1, p + 1)

  ## first calculate top left part for regression coefficients
  for (i in 1:p) {
    for (j in 1:p) {
      obsInfo[i, j] <- with(object, sum((1 + y/theta)*fitted.values*x[, i]*x[, j]/
                                          (1 + fitted.values/theta)^2))
    }
  }

  ## information for dispersion parameter
  obsInfo[(p + 1), (p + 1)] <- with(object, -sum(trigamma(theta + y) - trigamma(theta) -
                                                   1/(fitted.values + theta) +
                                                   (theta + y)/(theta + fitted.values)^2 -
                                                   1/(fitted.values + theta) +
                                                   1/theta))
  ## covariance between regression coefficients and dispersion
  for (i in 1:p) {
    obsInfo[(p + 1), i] <- with(object, -sum(((y - fitted.values)*fitted.values /
                                                ((theta + fitted.values)^2)) * x[, i]))
    obsInfo[i, (p + 1)] <- obsInfo[(p + 1), i]
  }
  ## return variance covariance matrix
  structure(ginv(obsInfo), dimnames = rep(list(c(rownames(vcov(object)), "theta")), 2))
}

