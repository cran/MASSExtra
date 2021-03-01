### 2d kde new method
### 

#' @importFrom utils methods
NULL

floor_os   <- function(x, o = 0, s = 1) o + s*floor((x - o)/s)
ceiling_os <- function(x, o = 0, s = 1) o + s*ceiling((x - o)/s)

#' One-dimensional Kernel Density Estimate
#' 
#' A pure R implementation of an approximate one-dimensional KDE, similar
#' to \code{\link[stats]{density}} but using a different algorithm not involving
#' \code{\link[stats]{fft}}.  Two extra facilities are provided, namely
#' (a) the kernel may be given either as a character string to select one of a 
#' number of kernel functions provided, or a user defined R function, and (b) the
#' kde may be fitted beyond the prescribed limits for the result, and folded back
#' to emulate the effect of having known bounds for the distribution.
#'
#' @param x   A numeric vector for which the kde is required or (in methods)
#'            an object of class \code{"kde_1d"}
#' @param bw  The bandwidth function.  Note this must be an R function.
#' @param kernel The kernel function, specified either as a character string or as
#'               an R function. Partial matching of the character string is allowed.
#' @param n  Integer, the number of equally-spaced values in the abscissa of the kde
#' @param limits  numeric vector of length 2.  Prescribed x-range limits for the 
#'                x-range of the result.  May be infinite, but infinite values will be
#'                pruned back to an appropriate value as determined by the data.
#' @param cut The number of bandwidths beyond the range of the input x-values to use
#' @param na.rm  Logical value: should any missing values in x be silently removed?
#' @param adjust numeric value: a multiplier to be applied to the computed bandwidth.
#' @param fold  Logical value: should the kde be estimated beyond the prescribed limits
#'              for the result and 'folded back' to emulate the effect of having known
#'              range boundaries for the underlying distribution?
#' @param las,col,xlab,ylab base graphics parameters
#' @param ... currently ignored, except in method functions
#'
#' @return A list of results specifying the result of the kde computation, of class \code{"kde_1d"}
#' @export
#'
#' @examples
#' set.seed(1234)
#' u <- runif(5000)
#' kdeu0 <- kde_1d(u, limits = c(-Inf, Inf))
#' kdeu1 <- kde_1d(u, limits = 0:1, kernel = "epan", fold = TRUE)
#' plot(kdeu0, col = 4)
#' lines(kdeu1, col = "dark green")
#' fun <- function(x) (0 < x & x < 1) + 0
#' curve(fun, add=TRUE, col = "grey", n = 1000)
kde_1d <- function(x, bw = bw.nrd0,
                   kernel = c("gaussian", 
                              "biweight", "cosine", "epanechnikov", "logistic", 
                              "optCosine", "rectangular", "squaredCosine",
                              "triangular", "tricube", "triweight", "uniform"),
                   n = 512, limits = c(rx[1] - cut*bw, rx[2] + cut*bw), 
                   cut = 3, na.rm = FALSE, adjust = 1, fold = FALSE, ...) {
  stopifnot(non_numeric_x = is.numeric(x),
            short_x = length(unique(x)) > 2)
  data_name <- deparse(substitute(x))
  if(is.function(bw)) {
    bw <- bw(x) * adjust
  }
  stopifnot(non_numeric_bw = is.numeric(bw),
            bad_bw = length(bw) == 1 && bw > 0)
  if(is.character(kernel)) {
    kernel <- match.arg(kernel)
    substring(kernel, 0, 1) <- toupper(substring(kernel, 0, 1))
    kernel <- get(paste0("kernel", kernel), mode = "function")
  } else {
    stopifnot(kernel_not_a_function = is.function(kernel),
              bad_arguments = identical(args(kernel), args(kernelGaussian)))
  }
  stopifnot(bad_n = is.numeric(n) && length(n) == 1 && (n %% 1 == 0) && n > 2)
  rx <- range(x, na.rm = na.rm)
  if(fold) {
    limits[1] <- max(limits[1], rx[1] - cut*bw)
    limits[2] <- min(limits[2], rx[2] + cut*bw)
    stopifnot(bad_limits = all(x >= limits[1]) && all(x <= limits[2]))
    dx <- diff(limits)/(n-1)
    lo <-   floor_os(limits[1] - min(diff(limits), 5*bw), limits[1], dx)
    up <- ceiling_os(limits[2] + min(diff(limits), 5*bw), limits[2], dx)
    n_star <- 1 + round((up - lo)/dx)
    xo <- seq(lo, up, by = dx)
    fr <- tabulate(1 + pmin(floor((x - lo)/dx), n_star - 1), nbins = n_star)
  } else {
    lo <- max(limits[1], rx[1] - 5*bw)
    up <- min(limits[2], rx[2] + 5*bw)
    stopifnot(bad_limits = all(x >= lo) && all(x <= up))
    dx <- (up - lo)/(n-1)
    xo <- seq(lo, up, by = dx)
    fr <- tabulate(1 + pmin(floor((x - lo)/dx), n - 1), nbins = n)
  }
  yo <- outer(xo, xo, function(x, y) kernel(x, y, bw)) %*% fr
  if(fold) {
    if(any(low <- (xo < limits[1]+dx/10))) {
      k <- which(low)
      kl <- k[length(k)]
      k0 <- 2*kl - k
      yo[k0] <- yo[k0] + yo[k] 
    }
    if(any(upp <- (limits[2]-dx/10 < xo))) {
      k <- which(upp)
      ku <- k[1]
      k0 <- 2*ku - k
      yo[k0] <- yo[k0] + yo[k] 
    }  
    keep <- kl:ku
    yo <- yo[keep]
    xo <- xo[keep]
    lo <- limits[1]
    up <- limits[2]
  }
  structure(list(x = xo, y = yo/(sum(yo)*dx), bw = bw, n = length(x),
                 lower = lo, upper = up, data_name = data_name),
            class = "kde_1d")
}

#' @rdname kde_1d
#' @export
print.kde_1d <- function(x, ...) {
  cat("A Kernel Density Estimate of class ", class(x), "\n\n")
  labs <- c("Response", "Sample size", "Range", "Bandwidth", "Resolution")
  values <- with(x, c(data_name, n,
                      paste0(signif(lower, 4), " < x < ",
                             signif(upper, 4)),
                      signif(bw, 4),
                      length(x)))
  cat(paste0(format(labs), ": ", values), sep = "\n")
  generics <- attr(utils::methods(class = "kde_1d"), "info")$generic
  cat("\nMethods exist for generics:",
      paste(generics, collapse = ", "), "\n")
  invisible(x)
}

#' @rdname kde_1d
#' @export
plot.kde_1d <- function(x, ..., col = "steel blue", las = 1,
                        xlab = bquote(x == italic(.(x$data_name))),
                        ylab = expression(kde(italic(x)))) {
  with(x, {
    plot(x = c(lower, x, upper), y = c(0, y, 0),
         type = "l", xlab = xlab, ylab = ylab, axes = FALSE,
         col = col, panel.first = grid(), las = las, ...)
    axis(1, las = las)
    axis(2, las = las)
  })
  invisible(x)
}

#' @import grDevices
NULL

#' A Two-dimensional Kernel Density Estimate
#' 
#' A pure R implementation of an approximate two-dimensional kde computation, where
#' the approximation depends on the x- and y-resolution being fine, i.e. the number
#' of both x- and y-points should be reasonably large, at least 256.  The coding
#' follows the same idea as used in \code{\link[MASS]{kde2d}}, but scales much better
#' for large data sets. 
#'
#' @param x,y Numeric vectors of the same length specified in any way acceptable
#'            to \code{\link[grDevices]{xy.coords}}.  In methods, \code{x} will be
#'            an object of class \code{"kde_2d"}
#' @param bw bandwidths. May be a numeric vector of length 1 or 2, or a function, 
#'           or list of two bandwidth computation functions.  Short entities will
#'           be repeated to length 1.  The first relates to the x-coordinate and
#'           the second to the y.
#' @param kernel As for \code{\link{kde_1d}} though 1 or 2 values may be specified
#'               relating to x- and y-coordinates respectively.  Short entities will
#'               be repeated to length 2
#' @param n positive integer vector of length 1 or 2 specifying the resolution required
#'          in the x- and y-coordinates respectively.  Short values will be repeated to
#'          length 2.
#' @param x_limits,y_limits Numeric vectors specifying the limits required for the result
#' @param cut The number of bandwidths beyond the x- and y-range limits for the resuls.
#' @param na.rm  Should missing values be silently removed?
#' @param adjust A factor to adjust both bandwidths to regulate smoothness
#' @param las,col,xlab,ylab base graphics parameters
#' @param ... currently ignored, except in method functions
#'
#' @return A list of results of class \code{"kde_2d"}.  The result may be used directly
#'         in \code{\link[graphics]{image}} or \code{\link[graphics]{contour}}.
#' @export
#'
#' @examples
#' krc <- with(Boston, {
#'   criminality <- log(crim)
#'   spaciousness <- sqrt(rm)
#'   kde_2d(criminality, spaciousness, n = 128, kernel = "biweight")
#' })
#' plot(krc, xlab = expression(italic(Criminality)), ylab = expression(italic(Spaciousness)))
#' contour(krc, add = TRUE, lwd = 0.5)
#' points(sqrt(rm) ~ log(crim), Boston, pch = ".", cex = 2, col = "dark green")
#' 
#' with(krc, persp(x, 10*y, 3*z, border="transparent", col = "powder blue",
#'                 theta = 30, phi = 15, r = 20, scale = FALSE, shade = TRUE, 
#'                 xlab = "Criminality", ylab = "Spaciousness", zlab = "density"))
kde_2d <- function(x, y = NULL, bw = list(x = bw.nrd0, y = bw.nrd0), 
                   kernel = c("gaussian", 
                              "biweight", "cosine", "epanechnikov", "logistic", 
                              "optCosine", "rectangular", "squaredCosine",
                              "triangular", "tricube", "triweight", "uniform"),
                   n = 256,
                   x_limits = c(rx[1] - cut*bw["x"], rx[2] + cut*bw["x"]), 
                   y_limits = c(ry[1] - cut*bw["y"], ry[2] + cut*bw["y"]),
                   cut = 1, na.rm = FALSE, adjust = 53/45, ...) {
  xy <- grDevices::xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  if(na.rm) {
    if(any(nax <- is.na(x)) | any(nay <- is.na(y))) {
      x <- x[!(nax | nay)]
      y <- y[!(nax | nay)]
    }
  }
  stopifnot(x_has_nas = !any(is.na(x)), 
            y_has_nas = !any(is.na(y)), 
            x_y_lengths_unequal = length(x) == length(y),
            bad_n = !is.integer(n) || all(n > 2))
  knames <- if(is.null(names(kernel))) c("x", "y") else names(kernel)
  if(is.character(kernel)) {
    kernel <- match.arg(kernel, several.ok = TRUE)
    kernel <- rep(kernel, length.out = 2)
    substring(kernel, 0, 1) <- toupper(substring(kernel, 0, 1))
    kernel <- paste0("kernel", kernel)
    # kernel <- lapply(kernel, get, mode = "function")
    kernel = list(x = get(kernel[1], mode = "function"),
                  y = get(kernel[2], mode = "function"))
  } else if(is.function(kernel)) {
    if(length(kernel) == 1) {
      kernel <- rep(list(kernel), 2)
    }
  }
  kernel <- setNames(kernel, knames)
  if(is.function(bw)) {
    bw <- rep(list(bw), 2)
  }
  if(is.null(names(bw))) {
    bw <- setNames(bw, c("x", "y"))
  }
  bw <- mapply(FUN = function(x, y) x(y), bw, list(x = x, y = y))*adjust
  nnames <- if(is.null(names(n))) c("x", "y") else names(n)
  n <- setNames(rep(n, length.out = 2), nnames)
  rx <- range(x)
  ry <- range(y)
  x_limits <- c(max(x_limits[1], rx[1]-5*bw["x"]), min(x_limits[2], rx[2]+5*bw["x"]))
  y_limits <- c(max(y_limits[1], ry[1]-5*bw["y"]), min(y_limits[2], ry[2]+5*bw["y"]))
  xo <- seq(x_limits[1], x_limits[2], length.out = n["x"])
  yo <- seq(y_limits[1], y_limits[2], length.out = n["y"])
  dx <- diff(x_limits)/(n["x"] - 1)
  dy <- diff(y_limits)/(n["y"] - 1)
  nxy <- prod(n)
  Fr <- tabulate(pmax(1, pmin(1 + floor((x - x_limits[1])/dx) +
                                n["x"]*floor((y - y_limits[1])/dy), nxy)),
                 nbins = nxy)
  dim(Fr) <- n
  Hx <- outer(xo, xo, kernel[["x"]], sd = bw["x"])
  Hy <- outer(yo, yo, kernel[["y"]], sd = bw["y"])
  z <- tcrossprod(Hx %*% Fr, Hy)
  z <- z/(sum(z)*dx*dy)
  structure(list(x = xo, y = yo, z = z, bw = bw, n_sample = length(x),
                 lower = c(x = min(xo), y = min(yo)), 
                 upper = c(x = max(xo), y = max(yo)),
                 n_out = n, data_name = with(xy, c(x = ifelse(is.null(xlab), "x", xlab),
                                                   y = ifelse(is.null(ylab), "y", ylab)))),
            class = "kde_2d")
}

#' @rdname kde_2d
#' @export
print.kde_2d <- function(x, ...) {
  cat("A Two-dimensional Kernel Density Estimate\n\n")
  labs <- c("Responses", "Sample size", "Ranges", "Bandwidths", "Resolution")
  values <- with(x, c(paste(data_name, collapse = ", "),
                      n_sample,
                      paste(paste0(signif(lower["x"], 4), " < x < ",
                                   signif(upper["x"], 4)),
                            paste0(signif(lower["y"], 4), " < y < ",
                                   signif(upper["y"], 4)),
                            sep = ", "),
                      paste(signif(bw, 4), collapse = ", "),
                      paste(length(x), length(y), sep = ", ")))
  cat(paste0(format(labs), ": ", values), sep = "\n")
  generics <- attr(methods(class = "kde_2d"), "info")$generic
  cat("\nMethods exist for generics:",
      paste(generics, collapse = ", "), "\n")
  invisible(x)
}

#' @rdname kde_2d
#' @export
plot.kde_2d <- function(x, ..., las = 1,
                        xlab = bquote(italic(.(x$data_name[["x"]]))),
                        ylab = bquote(italic(.(x$data_name[["y"]]))),
                        col = hcl.colors(50, "YlOrRd", rev = TRUE)) {
  with(x, {
    image(x, y, z, las = las, xlab = xlab, ylab = ylab, col = col, ...)
  })
  invisible(x)
}

#' @import demoKde
NULL

#' #' @rdname kde_1d
#' #' @export
#' kernelBiweight <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(7)*sd
#'   ifelse((z <- abs(x-mean)) < h, 15/16*(1 - (z/h)^2)^2/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelCosine <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(1/(1-8/pi^2))*sd
#'   ifelse((z <- abs(x-mean)) < h, pi/4*cos((pi*z)/(2*h))/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelEpanechnikov <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(5)*sd
#'   ifelse((z <- abs(x-mean)) < h, 3/4*(1 - (z/h)^2)/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelGaussian <- function(x, mean = 0, sd = 1)
#'   dnorm(x, mean = mean, sd = sd)
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelLogistic <- function(x, mean = 0, sd = 1)
#'   stats::dlogis(x, mean, sqrt(3)/pi*sd)
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelOptCosine <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(1/(1-8/pi^2))*sd
#'   ifelse((z <- abs(x-mean)) < h, pi/4*cos((pi*z)/(2*h))/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelRectangular <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(3)*sd
#'   ifelse(abs(x-mean) < h, 1/(2*h), 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelSquaredCosine <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(3/(1-6/pi^2))*sd
#'   ifelse((z <- abs(x-mean)) < h, cos(pi*z/(2*h))^2/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelTriangular <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(24)*sd/2
#'   ifelse((z <- abs(x-mean)) < h, (1 - z/h)/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelTricube <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(243/35)*sd
#'   ifelse((z <- abs(x - mean)) < h, 70/81*(1 - (z/h)^3)^3/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelTriweight <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(9)*sd
#'   ifelse((z <- abs(x-mean)) < h, 35/32*(1 - (z/h)^2)^3/h, 0)
#' }
#' 
#' #' @rdname kde_1d
#' #' @export
#' kernelUniform <- function(x, mean = 0, sd = 1) {
#'   h <- sqrt(3)*sd
#'   ifelse(abs(x-mean) < h, 1/(2*h), 0)
#' }
