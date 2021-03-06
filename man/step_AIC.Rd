% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stepAIC_BIC.R
\name{step_AIC}
\alias{step_AIC}
\alias{step_BIC}
\alias{step_GIC}
\alias{drop_term}
\alias{add_term}
\title{Stepwise model construction and inspection}
\usage{
step_AIC(object, ..., trace = 0, k = 2)

step_BIC(object, ..., trace = 0, k = max(2, log(nobs(object))))

step_GIC(object, ..., trace = 0, k = (2 + log(nobs(object)))/2)

drop_term(
  object,
  ...,
  test = default_test(object),
  k,
  sorted = TRUE,
  decreasing = TRUE,
  delta = TRUE
)

add_term(
  object,
  ...,
  test = default_test(object),
  k,
  sorted = TRUE,
  decreasing = TRUE,
  delta = TRUE
)
}
\arguments{
\item{object}{as for \code{\link[MASS]{stepAIC}}}

\item{...}{additional arguments passed on to main function in \code{MASS}}

\item{trace, k}{as for \code{\link[MASS]{stepAIC}}}

\item{sorted, test}{as for \code{\link[MASS]{dropterm}} and \code{\link[MASS]{addterm}}}

\item{decreasing}{in \code{drop_term} should the rows be displayed in decreasing order,
that is best to worst terms, from that of \code{\link[MASS]{dropterm}}?}

\item{delta}{Should the criterion be displayed (FALSE) or the change in the
in the criterion relative to the present model (TRUE)?}
}
\value{
A fitted model object after stepwise refinement, or a data frame with
        extra class membership for single term functions.
}
\description{
Front-ends to \code{\link[MASS]{stepAIC}} and \code{\link[MASS]{dropterm}} with changed defaults.
\code{step_BIC} implements a stepwise selection with BIC as the criterion and
\code{step_GIC} uses an experimental criterion with a penalty midway between AIC and BIC: the
"Goldilocks" criterion.
}
\examples{
fm <- glm.nb(Days ~ .^3, quine)
drop_term(fm_aic <- step_AIC(fm))
drop_term(fm_bic <- step_BIC(fm))
}
