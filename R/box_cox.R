#' @import methods
#' @import graphics
#' @import stats
NULL

#' Find the box-cox transform exponent estimate
#'
#' Estimates the box-cox power transformation appropriate for
#' a linear model
#'
#' @param bc either a \code{"box_cox"} object, a formula,data pair, a linear model object or an xy-lixt
#' @param data a data frame or envinonment
#' @param ... additional parameters passed on to \code{box_cox}
#' @param span integer: how many steps on either side of the maximum to use for the quadratic interpolation to find the maximum
#'
#' @return numeric: the maximum likelihood estimate of the exponent
#' @export
#'
#' @examples
#' lambda(medv ~ ., Boston, span = 10)
lambda <- function(bc, ...) {
  UseMethod("lambda")
}

#' @rdname lambda
#' @export
lambda.formula <- function(bc, data = sys.parent(), ..., span = 5) {
  call <- match.call(expand.dots = FALSE)
  call[[1]] <- as.name("lm")
  call$... <- call$span <- NULL
  names(call) <- sub("bc", "formula", names(call))
  bc <- eval(call, sys.parent())
  lambda.box_cox(box_cox(bc, ...), span = span)
}

#' @rdname lambda
#' @export
lambda.lm <- function(bc, ..., span = 5) {
  lambda.box_cox(box_cox(bc, ...), ..., span = span)
}

#' @rdname lambda
#' @export
lambda.box_cox <- function(bc, ..., span = 5) {
  with(bc, {
    n <- length(x)
    k <- which(y == max(y))[1]
    if(k < 1 + span || k > n - span)
      span <- 2*span
    is <- max(1, k - span):min(n, k + span)
    x0 <- x[is] - x[k]
    y0 <- y[is]
    b <- qr.coef(qr(cbind(1, x0, x0^2)), y0)
    if(b[3] >= 0)
      stop("the likelihood maximum is not at a stationary point")
    structure(as.vector(x[k] - b[2]/(2*b[3])),
              coefs = b, x0 = x[k], class = "lambda")
  })
}

#' @rdname lambda
#' @export
lambda.default <- function(bc, ...) {
  lambda.box_cox(bc, ...)
}

#' Print method for Box-Cox objects
#'
#' @param x an object of class \code{"box_cox"}
#' @param ... ignored
#'
#' @return x, invisibly
#' @export
print.lambda <- function(x, ...) {
  y <- x
  x <- as.vector(x)
  NextMethod()
  invisible(y)
}

## Lambda <- Vectorize(lambda, vectorize.args = "span")

#' Box-cox constructor function
#'
#' A front-end to \code{\link[MASS]{boxcox}} with slicker display and better defaults
#'
#' @param object either a \code{"box_cox"} object, a formula,data pair, a linear model object or an xy-lixt
#' @param data a data frame or environment
#' @param ... additional arguments passed on to methods
#' @param plotit currently ignored.  Plotting is done by \code{plot} or \code{print} methods
#' @param flap fraction of the central 95\% notional confidence to expand the range of lambda for the display
#' @param col.lines colour to use for indicator lines in the display
#' @param x a \code{"box_cox"} object to be displayed
#' @param xlab,ylab,las as for \code{plot}
#'
#' @return an object of class \code{"box_cox"}
#' @export
#'
#' @examples
#' box_cox(MPG.city ~ Weight, Cars93)
setGeneric("box_cox", function(object, ...) {
  standardGeneric("box_cox")
})

#' @rdname box_cox
#' @export
setMethod("box_cox", signature(object = "formula"),
          function(object, data = sys.parent(), ...) {
            Call <- substitute(lm(FORM, DATA),
                               list(FORM = object, DATA = substitute(data), qr = TRUE))
            object <- eval(Call, sys.parent(2))
            callGeneric(object, ...)
          })

#' @rdname box_cox
#' @export
setMethod("box_cox", signature(object = "lm"),
          function(object, ..., plotit, flap = 0.4) {
            bc <- boxcox(object, plotit = FALSE, ...)
            lam <- lambda(bc)
            b <- attr(lam, "coefs")
            x0 <- attr(lam, "x0")
            z <- qchisq(0.95, 1)/2
            lims <- sort(Re(polyroot(c(z, b[-1])))) + x0
            gap <- flap*diff(lims)
            out <- boxcox(object, plotit = FALSE,
                                lambda = seq(lims[1] - gap,
                                             lims[2] + gap,
                                             len = 500))
            out$y <- out$y - max(out$y)
            attr(out, "hlim") <- with(out, max(y)) - z
            attr(out, "lims") <- with(out, range(x[y >= max(y) - z]))
            attr(out, "lambda") <- lambda.default(out)
            class(out) <- "box_cox"
            out
          })

#' @rdname box_cox
#' @export
plot.box_cox <- function(x, ..., las = 1,
                         xlab = expression(lambda),
                         ylab, col.lines = "steel blue") {
  if(missing(ylab))
    ylab <- expression({X[LR]^2}(lambda) ==
                         2*`{`*log*' '*italic(L)(hat(lambda)) -
                         log*' '*italic(L)(lambda)*`}`)
  lims <- attr(x, "lims")
  hlim <- -2*attr(x, "hlim")
  lambda <- attr(x, "lambda")
  plot(-2*y ~ x, data.frame(unclass(x)), axes = FALSE,
       ann = FALSE, type = "l")
  pu <- par("usr")
  xx <- (pu[1] + lims[1])/2
  arrows(c(xx, xx), c(1.2, hlim-1.3), c(xx, xx), c(0, hlim),
         col = col.lines, length = 0.125, angle = 15)
  text(xx, hlim/2, expression(chi["0.95,1"]^2))
  abline(h = c(0, hlim), lty = "solid", col = col.lines)
  xx <- c(lims, lambda)
  yy <- c(hlim, hlim, max(x$y))
  arrows(xx, yy, xx, pu[4], ## , length = 0.125, angle = 15
         lty = "solid", col = col.lines, length = 0.125, angle = 15)
  axis(1)
  axis(2, las = las)
  at <- c(lims[1], lambda, lims[2])
  axis(3, at = at, labels = format(round(at, 3), digits = 3),
       cex.lab = 0.5, col.axis = col.lines)
  title(xlab = xlab, ylab = ylab, line = 2)
  box()
  with(unclass(x), lines(x, -2*y))
  grid()
  invisible(x)
}

#' @rdname box_cox
#' @export
print.box_cox <- plot.box_cox


#' Box-Cox transform
#'
#' Compute the box-cox transform of a vector of values, handling
#' the region near lambda = 0 with some care
#'
#' @param y numeric, the original observations
#' @param lambda numeric, the box-cox power
#' @param eps numeric, a guard aroung lambda = 0
#'
#' @return A vector of transformed quantities
#' @export
#'
#' @examples
#' plot(12:50, bc(12:50, -1), type = "l", xlab = "MPG", ylab = "bc(MPG, -1)",
#'      las = 1, col = "sky blue", panel.first = grid())
#' points(bc(MPG.city, -1) ~ MPG.city, data = Cars93, pch = 16, cex = 0.7)
bc <- function(y, lambda, eps = 1.0e-4) { ## vectorized wrt y and lambda
  n <- max(length(y), length(lambda))
  lambda <- rep_len(lambda, length.out = n)
  ly <- log(y)
  aly <- lambda*ly
  ifelse(abs(lambda) > eps,  (exp(aly) - 1)/lambda,
         ly*(1 + aly*(1 + aly*(1 + aly*(1 + aly*(1 + aly*(1 + aly/7)/6)/5)/4)/3)/2))
}

#' Box-Cox transform inverse
#'
#' Find the original value corresponding to a box-cox transform
#'
#' @param z numeric, the transformed value
#' @param lambda numeric, the power of the box-cox transform
#' @param eps numeric, a guard around lambda = 0
#'
#' @return A vector of original quantities
#' @export
#'
#' @examples
#' invy <- with(Cars93, bc(MPG.city, lambda = -1))
#' mpgc <- bc_inv(invy, lambda = -1)
#' range(mpgc - Cars93$MPG.city)
bc_inv <- function(z, lambda, eps = 1e-5) {
  n <- max(length(z), length(lambda))
  lambda <- rep_len(lambda, length.out = n)
  ifelse(abs(lambda) > eps, (1 + lambda*z)^(1/lambda),
         exp(z) * (1 - (z^2*lambda)/2 + ((3*z^4 + 8*z^3)*lambda^2)/24 -
                     ((z^6 + 8*z^5 + 12*z^4)*lambda^3)/48 +
                     ((15*z^8 + 240*z^7 + 1040*z^6 + 1152*z^5)*lambda^4)/5760 -
                     ((3*z^10 + 80*z^9 + 680*z^8 + 2112*z^7 + 1920*z^6)*lambda^5)/11520))
}
