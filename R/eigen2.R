### eigen2
###

#' Generalized eigenvalue problem
#'
#' Solves the generalized eigenvalue problem (B - lambda*W)*alpha = 0,
#' where B and W are symmetric matrices of the same size, W is positive
#' definite, lambda is a scalar and alpha and 0 are vectors.
#'
#' If W is not specified, W = I is assumed.
#'
#' @param B,W  Similarly sized symmetric matrices with W positive definite.
#'
#' @return A list with components \code{values} and \code{vectors} as for \code{\link[base]{eigen}}
#' @export
#'
#' @examples
#' X <- as.matrix(subset(iris, select = -Species))
#' W <- crossprod(resid(aov(X ~ Species, iris)))
#' B <- crossprod(resid(aov(X ~ 1,       iris))) - W
#' n <- nrow(iris)
#' p <- length(levels(iris$Species))
#' (ev <- eigen2(B/(p - 1), W/(n - p)))  ## hand-made discriminant analysis
#' DF <- X %*% ev$vectors[, 1:2]
#' with(iris, {
#'      plot(DF, col = Species, pch = 20,
#'           xlab = expression(DF[1]), ylab = expression(DF[2]))
#'      legend("topleft", levels(Species), pch = 20, col = 1:3)
#' })
eigen2 <- function(B, W) {
  if(!(is.numeric(B) && is.matrix(B) && nrow(B) == ncol(B)))
    stop(sQuote(deparse(substitute(B))), " is not a square numeric matrix.")

  if(missing(W)) return(eigen(B))

  if(!(is.numeric(W) && is.matrix(W) && nrow(W) == ncol(W)
       && isTRUE(all.equal(W, t(W)))))
    stop(sQuote(deparse(substitute(W))), " is not a numeric, symmetric matrix.")

  V <- try(chol(W), silent = TRUE)
  if(inherits(V, "try-error"))
    stop("The matrix ", sQuote(deparse(substitute(W))), " is not positive definite")
  V <- backsolve(V, diag(nrow(W)))
  tB <- crossprod(V, B %*% V)
  within.list(eigen(tB), vectors <- V %*% vectors)
}
