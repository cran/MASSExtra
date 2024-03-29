% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/avoid.R
\name{usr2in}
\alias{usr2in}
\alias{usr2in,numeric-method}
\alias{usr2in,xy-method}
\alias{in2usr}
\alias{in2usr,numeric-method}
\alias{in2usr,xy-method}
\title{Conversion functions for plotting}
\usage{
usr2in(x, ...)

\S4method{usr2in}{numeric}(
  x,
  y,
  usr = par("usr"),
  pin = par("pin"),
  xlog = par("xlog"),
  ylog = par("ylog"),
  ...
)

\S4method{usr2in}{xy}(x, ...)

in2usr(x, ...)

\S4method{in2usr}{numeric}(
  x,
  y,
  usr = par("usr"),
  pin = par("pin"),
  xlog = par("xlog"),
  ylog = par("ylog"),
  ...
)

\S4method{in2usr}{xy}(x, ...)
}
\arguments{
\item{x, y}{any of the forms that the coordinates of a scatterplot may be specified}

\item{...}{additional arguments for methods}

\item{usr, pin}{graphics parameters \code{par("usr"), par("pin")} (or replacements)}

\item{xlog, ylog}{logicals: are the x- and/or y-scales logarithmic?}
}
\value{
a \code{complex} vector of converted coordinates
}
\description{
Convert user coordinates to inch-based cordinates for the open display,
and back again
}
