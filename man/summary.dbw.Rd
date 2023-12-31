% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.dbw.R
\name{summary.dbw}
\alias{summary.dbw}
\title{Summarize the distribution balancing weighting estimation results}
\usage{
\method{summary}{dbw}(object, ...)
}
\arguments{
\item{object}{an object of class “dbw”, usually, a result of a call to \code{\link{dbw}}.}

\item{...}{additional arguments to be passed to summary.}
}
\value{
\item{call}{the matched call.}
\item{est}{the point estimate of the parameter of interest.}
\item{coefficients}{a table including coefficients, standard errors,
z-values, and two-sided p-values.}
\item{effn}{the effective sample size for the parameter of interest
estimation.}
}
\description{
Prints a summary of a fitted \code{dbw} object.
}
\details{
Prints a summary of a \code{dbw} object, in a format similar to glm.
}
\examples{
# For examples see example(dbw)
}
\seealso{
\code{\link{dbw}}, \code{\link[base]{summary}}
}
\author{
Hiroto Katsumata
}
