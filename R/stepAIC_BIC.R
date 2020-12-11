### step AIC and BIC
###

#' Stepwise model construction and inspection
#'
#' Front-ends to \code{\link[MASS]{stepAIC}} and \code{\link[MASS]{dropterm}} with changed defaults.
#' \code{step_BIC} implements a stepwise selection with BIC as the criterion and
#' \code{step_GIC} uses an experimental criterion with a penalty midway between AIC and BIC: the
#' "Goldilocks" criterion.
#'
#' @param object as for \code{\link[MASS]{stepAIC}}
#' @param ... additional arguments passed on to main function in \code{MASS}
#' @param trace,k as for \code{\link[MASS]{stepAIC}}
#' @param sorted,test as for \code{\link[MASS]{dropterm}} and \code{\link[MASS]{addterm}}
#' @param decreasing in \code{drop_term} should the rows be displayed in decreasing order,
#'        that is best to worst terms, from that of \code{\link[MASS]{dropterm}}?
#' @param delta Should the criterion be displayed (FALSE) or the change in the
#'        in the criterion relative to the present model (TRUE)?
#'
#' @return A fitted model object after stepwise refinement, or a data frame with
#'         extra class membership for single term functions.
#' @export
#'
#' @examples
#' fm <- glm.nb(Days ~ .^3, quine)
#' drop_term(fm_aic <- step_AIC(fm))
#' drop_term(fm_bic <- step_BIC(fm))
step_AIC <- function(object, ..., trace = 0, k = 2) {
  if(isTRUE(object$call$trace)) {
    warning("Trace detected. Turning it off.", immediate. = TRUE)
    call <- substitute(update(OBJ, trace = FALSE),
                       list(OBJ = substitute(object)))
    object <- eval.parent(call)
  }
  out <- stepAIC(object, ..., trace = trace, k = k)
  attr(out, "penalty") <- k
  out
}

#' @rdname step_AIC
#' @export
step_BIC <- function(object, ..., trace = 0,
                    k = max(2, log(nobs(object)))) {
  if(isTRUE(object$call$trace)) {
    warning("Trace detected. Turning it off.", immediate. = TRUE)
    call <- substitute(update(OBJ, trace = FALSE),
                       list(OBJ = substitute(object)))
    object <- eval.parent(call)
  }
  out <- stepAIC(object, ..., trace = trace, k = k)
  attr(out, "penalty") <- k
  out
}

#' @rdname step_AIC
#' @export
step_GIC <- function(object, ..., trace = 0,
                     k = (2 + log(nobs(object)))/2) {
  if(isTRUE(object$call$trace)) {
    warning("Trace detected. Turning it off.", immediate. = TRUE)
    call <- substitute(update(OBJ, trace = FALSE),
                       list(OBJ = substitute(object)))
    object <- eval.parent(call)
  }
  out <- stepAIC(object, ..., trace = trace, k = k)
  attr(out, "penalty") <- k
  out
}

#' @rdname step_AIC
#' @export
drop_term <- function(object, ..., test = default_test(object), k, 
                      sorted = TRUE, decreasing = TRUE, delta  = TRUE) {
  if(!isS4(object) && isTRUE(object$call$trace)) {
    warning("Trace detected. Turning it off.", immediate. = TRUE)
    call <- substitute(update(OBJ, trace = FALSE),
                      list(OBJ = substitute(object)))
    object <- eval.parent(call)
  }
  k <- if(missing(k)) {
    p <- attr(object, "penalty")
    ifelse(is.null(p), 2, p)
  } else {
    if(is.character(k)) {
      switch(tolower(k),
             bic = max(2, log(nobs(object))),
             gic = (2 + log(nobs(object)))/2,
             aic =, 2)
    } else k
  }
  out <- MASS::dropterm(object, ..., test = test, k = k)
  if(sorted) out <- out[order(out$AIC, decreasing = decreasing), ]
  if(delta) {
    out$AIC <- out$AIC - out$AIC[rownames(out) == "<none>"]
    names(out) <- sub("AIC", "delta_AIC", names(out))
  }
  if(k != 2) {
    isBIC <- isTRUE(all.equal(k, log(nobs(object))))
    isGIC <- isTRUE(all.equal(k, (2 + log(nobs(object)))/2))
    names(out) <- sub("AIC", ifelse(isBIC, "BIC",
                                    ifelse(isGIC, "GIC",
                                           paste0("IC(", round(k, 2), ")"))),
                      names(out))
  }
  structure(out, pName = sub("^delta_", "", grep("IC", names(out), value = TRUE)[1]),
            class = c("drop_term", class(out)))
}

#' @rdname step_AIC
#' @export
add_term <- function(object, ..., test = default_test(object), k, 
                      sorted = TRUE, decreasing = TRUE, delta  = TRUE) {
  if(!isS4(object) && isTRUE(object$call$trace)) {
    warning("Trace detected. Turning it off.", immediate. = TRUE)
    call <- substitute(update(OBJ, trace = FALSE),
                       list(OBJ = substitute(object)))
    object <- eval.parent(call)
  }
  k <- if(missing(k)) {
    p <- attr(object, "penalty")
    ifelse(is.null(p), 2, p)
  } else {
    if(is.character(k)) {
      switch(tolower(k),
             bic = max(2, log(nobs(object))),
             gic = (2 + log(nobs(object)))/2,
             aic =, 2)
    } else k
  }
  out <- MASS::addterm(object, ..., test = test, k = k)
  if(sorted) out <- out[order(out$AIC, decreasing = decreasing), ]
  if(delta) {
    out$AIC <- out$AIC - out$AIC[rownames(out) == "<none>"]
    names(out) <- sub("AIC", "delta_AIC", names(out))
  }
  if(k != 2) {
    isBIC <- isTRUE(all.equal(k, log(nobs(object))))
    isGIC <- isTRUE(all.equal(k, (2 + log(nobs(object)))/2))
    names(out) <- sub("AIC", ifelse(isBIC, "BIC",
                                    ifelse(isGIC, "GIC",
                                           paste0("IC(", round(k, 2), ")"))),
                      names(out))
  }
  structure(out, pName = sub("^delta_", "", grep("IC", names(out), value = TRUE)[1]),
            class = c("drop_term", "add_term", class(out)))
}

#' drop_term plot method
#'
#' @param x An object of class \code{"drop_term"} generated by tither \code{drop_term} or \code{add_term}
#' @param ...,horiz arguments past on to \code{graphics::barplot}
#' @param col,border \code{barplot} fill and border colour(s) for positive and negative changes to the criterion, respectively
#' @param las graphics parameter
#' @param show.model logical: should the model itself be displayed?
#'
#' @return \code{x} invisibly
#' @export
#' @examples
#' boston_quad <- lm(medv ~ . + (rm + tax + lstat)^2 + poly(rm, 2) + 
#'          poly(tax, 2) + poly(lstat, 2), Boston)
#' dboston_quad <- drop_term(boston_quad, k = "bic")
#' plot(dboston_quad)
#' plot(dboston_quad, horiz = FALSE)
plot.drop_term <- function(x, ..., horiz = TRUE,
                          las = ifelse(horiz, 1, 2),
                          col = c("#DF536B", "#2297E6"),
                          border = c("#DF536B", "#2297E6"),
                          show.model = TRUE) {
  pName <- attr(x, "pName")
  stopifnot(isTRUE(any(grepl(pName, names(x)))))
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  # names(x)[grep(pName, names(x))[1]] <- pName ## remove any decoration such as "delta"
  AIC <- x[[which(grepl(pName, names(x)))[1]]]
  names(AIC) <- rownames(x)
  AIC <- sort(AIC - AIC["<none>"])
  if(is.character(col) && length(col) == 2) {
    col <- ifelse(AIC < 0, col[1], col[2])
  }
  if(is.character(border) && length(border) == 2) {
    border <- ifelse(AIC < 0, border[1], border[2])
  }
  pmar <- pmax(par("mar"), (if(horiz) c(4,6,1,1) else c(6,4,1,1)) + 0.1)
  pmar[3:4] <- 1.1
  par(mar = pmar, cex.axis = 0.8)
  if(show.model) {
    h <- attr(x, "heading")
    h <- format(gsub("\n", "", h[!grepl("^Single ", h)]), justify = "left")
    layout(cbind(2:1), heights = c(1, 6), widths = 1, respect = FALSE)
    on.exit(layout(matrix(1)), add = TRUE)
  }
  if(horiz) {
    barplot(AIC, xlab = bquote(Delta*' '*.(pName)), horiz = TRUE,
            las = las, col = col, border = border, ...)
  } else {
    barplot(AIC, ylab = bquote(Delta*' '*.(pName)), horiz = FALSE,
            las = las, col = col, border = border, ...)
  }
  if(show.model) {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("bottomleft", legend = h, bty = "n", cex = 0.8)
  }
  invisible(x)
}

#' Guess the default test
#'
#' Find an appropriate test to use in \code{\link[MASS]{dropterm}} if not specified
#'
#' @param object a fitted model object accommodated by \code{\link[MASS]{dropterm}}
#'
#' @return A character string, one of \code{"F", "Chisq", or "none"}
#' @export
#'
#' @examples
#' fm <- glm.nb(Days ~ .^3, quine)
#' default_test(fm)
default_test <- function(object) {
  UseMethod("default_test")
}

#' @rdname default_test
#' @export
default_test.default <- function(object) {
  "none"
}

#' @rdname default_test
#' @export
default_test.negbin <- function(object) {
  "Chisq"
}

#' @rdname default_test
#' @export
default_test.lmerMod <- function(object) {
  "Chisq"
}

#' @rdname default_test
#' @export
default_test.glmerMod <- function(object) {
  "Chisq"
}

#' @rdname default_test
#' @export
default_test.multinom <- function(object) {
  ## same as 'negbin'
  "Chisq"
}

#' @rdname default_test
#' @export
default_test.polr <- function(object) {
  ## same as 'negbin'
  "Chisq"
}

#' @rdname default_test
#' @export
default_test.glm <- function(object) {
  switch(object$family$family,
         binomial =, poisson = "Chisq",
         gaussian =, Gamma =, quasi =, quasibinomial =,
         quasipoisson =, inverse.gaussian = "F",
         "none")
}

#' @rdname default_test
#' @export
default_test.lm <- function(object) {
  "F"
}


############################################################################
############################################################################
###                                                                      ###
###             A CRUDE BACKWARD ELIMINATION TOOL FOR MODELS             ###
###                (GOOD ENOUGH FOR GOVERNMENT PURPOSES)                 ###
###                                                                      ###
############################################################################
############################################################################


.eliminate <- function(object, ..., trace) {
  d <- MASS::dropterm(object, sorted = TRUE, ...)
  if(trace)
    print(d)
  if(rownames(d)[1] %in% c("<none>", "1")) { ## hit the limit, do nothing
    return(invisible(object))
  }
  var <- parse(text = rownames(d)[1])[[1]]
  if(trace)
    cat("\nEliminating: ", rownames(d)[1], "\n-----------\n")
  upd <- substitute(update(OBJ, . ~ . - VAR), list(OBJ = substitute(object),
                                                   VAR = var))
  return(invisible(eval.parent(upd)))
}

#' Naive backeward elimination
#'
#' A simple facility to refine models by backward elimination.
#' Covers cases where \code{\link{drop_term}} works but \code{\link{step_AIC}}
#' does not
#'
#' @param object A fitted model object
#' @param ... additional arguments passed to \code{\link{drop_term}} such as \code{k}
#' @param trace logical: do you want a trace of the process printed?
#' @param k penalty (default 2, as for AIC)
#'
#' @return A refined fitted model object
#' @export
#'
#' @examples
#' fm <- lm(medv ~ . + (rm + tax + lstat)^2 +
#'            I((rm - 6)^2) + I((tax - 400)^2) + I((lstat - 12)^2), Boston)
#' sfm <- step_down(fm, trace = TRUE, k = "bic")
step_down <- function(object, ..., trace = FALSE, k) {
  if(!isS4(object) && isTRUE(object$call$trace)) {
    warning("Trace detected. Turning it off.", immediate. = TRUE)
    call <- substitute(update(OBJ, trace = FALSE),
                       list(OBJ = substitute(object)))
    object <- eval.parent(call)
  }
  k <- if(missing(k)) {
    p <- attr(object, "penalty")
    ifelse(is.null(p), 2, p)
  } else {
    if(is.character(k)) {
      switch(tolower(k),
             bic = max(2, log(nobs(object))),
             gic = (2 + log(nobs(object)))/2,
             aic =, 2)
    } else k
  }
  # oldOpt <- options(warn = -1)
  # on.exit(options(oldOpt))

  obj <- object
  suppressWarnings({
    repeat {
      tmp <- .eliminate(obj, ..., trace = trace, k = k)
      if(identical(tmp, obj)) break
      obj <- tmp
    }
  })
  attr(obj, "penalty") <- k
  obj
}


#' Intermediate Information Criterion
#'
#' An AIC-variant criterion that weights complexity with a penalty
#' mid-way between 2 (as for AIC) and log(n) (as for BIC).  I.e.
#' "not too soft" and "not too hard", just "Glodilocks".
#'
#' @param object a fitted model object for which the criterion is desired
#'
#' @return The GIC criterion value
#' @export
#'
#' @examples
#' gm <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), quine)
#' c(AIC = AIC(gm), GIC = GIC(gm), BIC = BIC(gm))
GIC <- function(object) {
  stats::AIC(object, k = (2 + log(stats::nobs(object)))/2)
}
