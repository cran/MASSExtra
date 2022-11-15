#' An S4 class to represent alternavive complex, matrix or list input forms.
#' @export
setClassUnion("xy", c("complex", "matrix", "list"))

setValidity("xy", function(object) {
  if(is.complex(object)) return(TRUE)
  if(is.list(object)) {
    if(all(c("x", "y") %in% names(object)) &&
       is.numeric(object$x) &&
       is.numeric(object$y)) {
      return(TRUE)
    } else {
      return("Invalid list argument.")
    }
  }
  if(is.matrix(object)) {
    if(is.numeric(object) && ncol(object) == 2) {
      return(TRUE)
    } else {
      return("Invalid matrix object")
    }
  }
})


#' Coerce to complex
#'
#' Utility function to create complex vectors from arguments
#' specified as in grDevices::xy.coords() or otherwise
#'
#' @param x A numeric vector or missing, or an object inheriting from class "xy"
#' @param y If x is a numeric an optional numeric vector, or missing. If x or y are
#'    missing they are taken as 0, but only one may be missing.
#'
#' @return A complex vector specifying 2-dimensional coordinates
#' @export
#'
#' @examples
#' as_complex(cbind(1:3, 3:1))
#' as_complex(y = 1:3)  ## real parts all zero
setGeneric("as_complex", function(x, y) {
  standardGeneric("as_complex")
})

#' @rdname as_complex
#' @export
setMethod("as_complex", signature(x = "xy", y = "missing"),
          function(x) {
            if(isTRUE(check <- validObject(x, TRUE))) {
              with(grDevices::xy.coords(x, recycle = TRUE), complex(real = x, imaginary = y))
            } else {
              stop(check)
            }
          })

#' @rdname as_complex
#' @export
setMethod("as_complex", signature(x = "numeric", y = "numeric"),
          function(x, y) {
            complex(real = x, imaginary = y)
          })

#' @rdname as_complex
#' @export
setMethod("as_complex", signature(x = "numeric", y = "missing"),
          function(x, y) {
            complex(real = x, imaginary = 0)
          })

#' @rdname as_complex
#' @export
setMethod("as_complex", signature(x = "missing", y = "numeric"),
          function(x, y) {
            complex(real = 0, imaginary = y)
          })

#' Avoid overlaps
#'
#' Generate a vector of positions to use to minimise text overlaps in labelled scatterplots
#'
#' @param x,y any of the forms that the coordinates of a scatterplot may be specified
#' @param ... additional arguments for methods
#' @param xlog,ylog logicals: are the x- and/or y-scales logarithmic?
#' @param usr,pin graphics parameters \code{par("usr"), par("pin")} (or replacements)
#' @param eps numeric: a zero tolerance
#' @param pi numeric: the value of the arithmetic constant of the same name
#'
#' @return a vector of integers all of which are 1, 2, 3, or 4, indicating placement positions.
#' @export
#'
#' @examples
#' set.seed(123)
#' z <- complex(real = runif(50), imaginary = runif(50))
#' mz <- mean(z)
#' z <- z[order(Arg(z - mz))]
#' plot(z, axes = FALSE, ann = FALSE)
#' segments(Re(mz), Im(mz), Re(z), Im(z))
#' abline(h = Im(mz), v = Re(mz), lwd = 0.5)
#' box()
#' text(Re(z), Im(z), pos = avoid(z), cex = 0.7, offset = 0.25,
#'      col = "red", font = 2, xpd = NA)
setGeneric("avoid", function(x, ...)
           standardGeneric("avoid"))

#' @rdname avoid
#' @export
setMethod("avoid", signature(x = "numeric"),
          function (x, y, ..., xlog = par("xlog"), ylog = par("ylog"),
                    usr = par("usr"), pin = par("pin"),
                    eps = .Machine$double.eps, pi = base::pi) {
            n <- max(length(x), length(y))
            z <- complex(real =      rep_len(x, length.out = n),
                         imaginary = rep_len(y, length.out = n))
            z <- usr2in(z, usr = usr, pin = pin)
            xydist <- outer(z, z, function(x, y) Mod(x - y))
            diag(xydist) <- Inf
            nearby <- apply(xydist, 2, which.min)
            zdiff <- z - z[nearby]
            pos <- findInterval(-(Arg(zdiff) + pi/4) %% (2*pi),
                                pi/2*(0:4), all.inside = TRUE)
            for(k in which(Mod(zdiff) <= eps)) {
              pos[sort(c(k, nearby[k]))] <- c(3, 1)
            }
            setNames(pos, seq_along(pos))
          })

#' @rdname avoid
#' @export
setMethod("avoid", signature(x = "xy"),
          function(x, ...) {
            xy <- grDevices::xy.coords(x, recycle = TRUE)
            out <- callGeneric(xy$x, xy$y, ...)
            if(length(names(x)) == length(out)) names(out) <- names(x)
            out
          })

#' Conversion functions for plotting
#'
#' Convert user coordinates to inch-based cordinates for the open display,
#' and back again
#'
#' @param x,y any of the forms that the coordinates of a scatterplot may be specified
#' @param ... additional arguments for methods
#' @param xlog,ylog logicals: are the x- and/or y-scales logarithmic?
#' @param usr,pin graphics parameters \code{par("usr"), par("pin")} (or replacements)
#'
#' @return a \code{complex} vector of converted coordinates
#' @export
setGeneric("usr2in", function(x, ...) {
  standardGeneric("usr2in")
})

#' @rdname usr2in
#' @export
setMethod("usr2in", signature(x = "numeric"),
          function(x, y, usr = par("usr"), pin = par("pin"),
                   xlog = par("xlog"), ylog = par("ylog"), ...) {
            xy <- grDevices::xy.coords(x, y, recycle = TRUE) ## safety
            with(xy, {
              x <- ((if(xlog) log(x) else x) - usr[1])/diff(usr[1:2])*pin[1]
              y <- ((if(ylog) log(y) else y) - usr[3])/diff(usr[3:4])*pin[2]
              complex(real = x, imaginary = y)
            })
          })

#' @rdname usr2in
#' @export
setMethod("usr2in", signature(x = "xy"),
          function(x, ...) {
            xy <- grDevices::xy.coords(x, recycle = TRUE)
            callGeneric(xy$x, xy$y, ...)
          })

#' @rdname usr2in
#' @export
setGeneric("in2usr", function(x, ...)
  standardGeneric("in2usr"))

#' @rdname usr2in
#' @export
setMethod("in2usr", signature(x = "numeric"),
          function(x, y, usr = par("usr"), pin = par("pin"),
                   xlog = par("xlog"), ylog = par("ylog"), ...) {
            xy <- grDevices::xy.coords(x, y, recycle = TRUE) ## safety
            with(xy, {
              x <- usr[1] + x/pin[1]*diff(usr[1:2])
              y <- usr[3] + y/pin[2]*diff(usr[3:4])
              complex(real =      if(xlog) exp(x) else x,
                      imaginary = if(ylog) exp(y) else y)
            })
          })

#' @rdname usr2in
#' @export
setMethod("in2usr", signature(x = "xy"),
          function(x, ...) {
            xy <- grDevices::xy.coords(x, recycle = TRUE)
            callGeneric(xy$x, xy$y, ...)
          })

#' Unit change functions
#'
#' Convert imperial to metric units, and vice versa.
#'
#' @param cm,inch,mm numeric vectors in the appropriate units
#'
#' @return a numeric vector of values in the new units
#' @name unitChange
NULL

#' @rdname unitChange
#' @export
cm2in <- function(cm) {
  cm / 2.54
}

#' @rdname unitChange
#' @export
mm2in <- function(mm) {
  mm / 25.4
}

#' @rdname unitChange
#' @export
in2cm <- function(inch) {
  inch * 2.54
}

#' @rdname unitChange
#' @export
in2mm <- function(inch) {
  inch * 25.4
}

