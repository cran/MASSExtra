Seq <- function(m, n) {
  unlist(mapply(FUN = seq, m, n))
}

#### lower triangle indices
####
#' Which in lower/upper triangle
#'
#' Find where the original positions of components are
#' in a matrix given a logical vector corresponding
#' to the lower or upper triangle stored by columns.  Similar
#' to which(.., arr.ind = TRUE)
#'
#' @param cond logical vector of length that of the lower triangle
#' @param diag logical: are the diagonal entries included?
#' @param lower logical: is this the lower triangle?  If FALSE it is the upper.
#'
#' @return a two column matrix with the row and column indices as the rows
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(20*2), 20, 2)
#' plot(X, asp = 1, pch = 16, las = 1, xlab = "x", ylab = "y")
#' dX <- dist(X)
#' ij <- which_tri(dX == max(dX))
#' points(X[as.vector(ij), ], col = "red", cex = 2, pch = 1)
#' segments(X[ij[1], 1], X[ij[1], 2],
#'          X[ij[2], 1], X[ij[2], 2], col = "red")
#' ij <- which_tri(dX == sort(dX, decreasing = TRUE)[2])
#' points(X[as.vector(ij), ], col = "blue", cex = 2, pch = 1)
#' segments(X[ij[1], 1], X[ij[1], 2],
#'          X[ij[2], 1], X[ij[2], 2], col = "blue")
#' polygon(X[chull(X), ], border = "sky blue")
#' rm(X, dX, ij)
which_tri <- function(cond, diag = FALSE, lower = TRUE) {
  stopifnot(is.logical(cond), length(cond) > 0)
  stopifnot(is.logical(diag), length(diag) == 1)
  stopifnot(is.logical(lower), length(lower) == 1)
  m <- length(cond)
  n <- (sqrt(8*m + 1) + 1)/2 - diag
  if(!isTRUE(all.equal(n, as.integer(n))) && (n > 1 - diag)) {
    stop("The length of cond is ", m,
         " which does not correspond to a lower triangle.")
  }
  if(lower) {
    cbind(row = Seq((2-diag):n,              n)[cond],
          col = rep(1:(n-1+diag), (n-1+diag):1)[cond])
  } else {
    cbind(row = Seq(1,          1:(n-1+diag))[cond],
          col = rep((2-diag):n, 1:(n-1+diag))[cond])
  }
}


eps <- .Machine$double.eps^0.75

rot <- function(x) {
  if(max(abs(x)) < eps) return(diag(2))
  if(abs(x[1]) > abs(x[2])) {
    t <- x[2]/x[1]
    a <- 1/sqrt(1 + t^2)
    b <- t*a
  } else {
    t <- x[1]/x[2]
    b <- 1/sqrt(1 + t^2)
    a <- t*b
  }
  rbind(c( a,  b),
        c(-b,  a))
}

#' Givens orthogonalisation
#'
#' Orthogonalization using Givens' method.
#'
#' @param X a numeric matrix with ncol(X) <= nrow(X)
#' @param nullspace logical: do you want an orthogonal basis for the null space?
#'
#' @return A list with components Q, R, as normally defined, and if nullspace is TRUE
#'    a further component N giving the basis for the requested null space of X
#' @export
#'
#' @examples
#' set.seed(1234)
#' X <- matrix(rnorm(7*6), 7)
#' givens_orth(X, nullspace = TRUE)
givens_orth <- function(X, nullspace = FALSE) {
  p <- ncol(X)
  n <- nrow(X)
  stopifnot(p <= n)
  A <- cbind(X, diag(n))
  j_to_P <- 1:ncol(A)
  for(j in 1:p) {
    i <- c(n-1, n)
    while(i[2] > j) {
      A[i,  j_to_P] <- rot(A[i, j]) %*% A[i, j_to_P]
      A[i[2], j] <- 0
      i <- i - 1
    }
    j_to_P <- j_to_P[-1]
  }
  if(any(neg <- diag(A) < 0)) {
    A[neg, ] <- -A[neg, ]
  }
  out <- list(Q = t(A[1:p, -(1:p), drop = FALSE]),
              R =   A[1:p,   1:p,  drop = FALSE])
  if(any(nullspace)) {
    if(p < n) {
      out$N <- t(A[-(1:p), -(1:p), drop = FALSE])
    } else {
      out$N <- matrix(0, n, 0)
    }
  }
  out
}

#' Gram-Schmidt orthogonalization
#'
#' Either classical or modified algorithms.
#' The modified algorithm is the more accurate.
#'
#' @param X a numerical matrix with ncol(X) <= nrow(X)
#'
#' @return A list with two components, Q, R, as usually defined.
#' @export
#'
#' @examples
#' set.seed(1234)
#' X <- matrix(rnorm(10*7), 10)
#' gs_orth_modified(X)
#' all.equal(gs_orth(X), gs_orth_modified(X))
#' all.equal(gs_orth_modified(X), givens_orth(X))
gs_orth_modified <- function(X) {
  n <- ncol(X)
  stopifnot(nrow(X) >= n)
  R <- matrix(0, n, n)
  for(j in 1:n) {
    v <- X[, j] <- X[,j]/(R[j,j] <- sqrt(sum(X[,j]^2)))
    if (j < n)
      for(i in (j+1):n) {
        X[, i] <- X[, i] - (R[j,i] <- sum(X[,i] * v)) * v
      }
  }
  list(Q = X, R = R)
}

#' @rdname gs_orth_modified
#' @export
gs_orth <- function(X) {
  n <- ncol(X)
  stopifnot(nrow(X) >= n)
  Q <- X
  R <- matrix(0, n, n)
  for(j in 1:n) {
    v <- X[, j]
    for(i in seq_len(j-1)) {
      v <- v - (R[i, j] <- sum(Q[, i]*X[, j]))*Q[, i]
    }
    Q[, j] <- v/(R[j,j] <- sqrt(sum(v^2)))
  }
  list(Q = Q, R = R)
}
