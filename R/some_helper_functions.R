
#' @import splines
#' @export
splines::ns
#' @export
splines::bs


#' Normalise a vector
#'
#' Similar to base::scale() but returning a vector with class attribute.
#' Used for safe prediction
#'
#' @param x A numeric vector
#' @param location A numeric vector of length 1
#' @param scale A numeric vector of length 1, usually positive
#'
#' @return A normalised vector inheriting from class "normalise"
#' @export
.normalise <- function(x, location, scale) {
  z <- structure(cbind(z = (x - location)/scale),
                 parms = c(location = unname(location),
                           scale = unname(scale)))
  class(z) <- c("normalise", "matrix")
  z
}

#' Method function for safe prediction
#'
#' This is an internal function not intended to be called directly by the user.
#'
#' @param var A numeric variable
#' @param call A single term from a linear model formula
#'
#' @return A call object used in safe prediction
#' @export
makepredictcall.normalise <- function(var, call) {
  parms <- attr(var, "parms")
  call <- call[1:2]
  call[[1]] <- as.name(".normalise")
  call[names(parms)] <- parms
  call
}

#' Standardisation functions for models
#'
#' These functions are for use in fitting linear models (or allies) with scaled
#' predictors, in such a way that when the fitted model objects are used for
#' prediction (or visualisation) the same scaling parameters will be used with
#' the new data.
#'
#' @param x A numeric vector
#'
#' @return a standardised vector containing the parameters needed for use in prediction with new data
#' @export
#'
#' @examples
#' fm <- lm(Gas ~ Insul/zs(Temp), whiteside)
#' gm <- lm(Gas ~ Insul/zu(Temp), whiteside)
#' hm <- lm(Gas ~ Insul/Temp,     whiteside)
#' c(fm = unname(predict(fm, data.frame(Insul = "Before", Temp = 0.0))),
#'   gm = unname(predict(gm, data.frame(Insul = "Before", Temp = 0.0))),
#'   hm = unname(predict(hm, data.frame(Insul = "Before", Temp = 0.0))))
#' rm(fm, gm, hm)
zs <- function(x) {  ## scale to z-scores, (equivalent to scale() in base)
  .normalise(x, mean(x), sd(x))
}

#' @rdname zs
#' @export
zu <- function(x) {  ## scale range to [0, 1]
  .normalise(x, min(x), diff(range(x)))
}

#' @rdname zs
#' @export
zr <- function(x) {  ## a 'robust' standardization
  .normalise(x, median(x), mad(x))
}

#' @rdname zs
#' @export
zq <- function(x) {  ## scale the box of a boxplot to [0, 1]
  .normalise(x, quantile(x, 0.25), IQR(x))
}

#' Mean and variance for a circular sample
#'
#' @param theta A vector of angles (in radians)
#'
#' @return The mean (rsp. variance) of the angle sample
#' @export
#'
#' @examples
#' th <- 2*base::pi*(rbeta(2000, 1.5, 1.5) - 0.5)
#' c(mn = mean_c(th), va = var_c(th))
#' rm(th)
mean_c <- function(theta) {
  Arg(mean(complex(argument = theta)))
}

#' @rdname mean_c
#' @export
var_c <- function(theta) {
  1 - Mod(mean(complex(argument = theta)))
}
