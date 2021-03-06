% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normal_sym_trunc_ab.R
\name{Normal_sym_trunc_ab}
\alias{Normal_sym_trunc_ab}
\alias{dNormal_sym_trunc_ab}
\alias{pNormal_sym_trunc_ab}
\alias{qNormal_sym_trunc_ab}
\alias{rNormal_sym_trunc_ab}
\alias{eNormal_sym_trunc_ab}
\alias{lNormal_sym_trunc_ab}
\title{The symmetric truncated normal distribution.}
\usage{
dNormal_sym_trunc_ab(
  x,
  sigma = 0.3,
  a = 0,
  b = 1,
  params = list(sigma, a, b),
  ...
)

pNormal_sym_trunc_ab(
  q,
  sigma = 0.3,
  a = 0,
  b = 1,
  params = list(mu = 2, sigma = 5, a = 0, b = 1),
  ...
)

qNormal_sym_trunc_ab(
  p,
  sigma = 0.3,
  a = 0,
  b = 1,
  params = list(mu = 2, sigma = 5, a = 0, b = 1),
  ...
)

rNormal_sym_trunc_ab(
  n,
  mu = 2,
  sigma = 3,
  a = 0,
  b = 1,
  params = list(sigma, a, b),
  ...
)

eNormal_sym_trunc_ab(X, w, method = "numerical.MLE", ...)

lNormal_sym_trunc_ab(
  X,
  w,
  mu = 2,
  sigma = 3,
  a = 0,
  b = 1,
  params = list(sigma, a, b),
  logL = TRUE,
  ...
)
}
\arguments{
\item{x, q}{A vector of quantiles.}

\item{a, b}{Boundary parameters.}

\item{params}{A list that includes all named parameters.}

\item{...}{Additional parameters}

\item{p}{A vector of probabilities.}

\item{n}{Number of observations.}

\item{mu, sigma}{Shape parameters.}

\item{X}{Sample observations.}

\item{w}{An optional vector of sample weights.}

\item{method}{Parameter estimation method.}

\item{logL}{logical;if TRUE, lNormal_sym_trunc_ab gives the log-likelihood, otherwise the likelihood is given.}
}
\value{
dNormal_sym_trunc_ab gives the density, pNormal_sym_trunc_ab the distribution function,
qNormal_sym_trunc_ab the quantile function, rNormal_sym_trunc_ab generates random deviates,and
eNormal_sym_trunc_ab estimates the parameters. lNormal_sym_trunc_ab provides the log-likelihood function.
}
\description{
Density, distribution, quantile, random number
generation and parameter estimation functions for the symmetric truncated normal distribution with parameters, \code{sigma},
\code{a} and \code{b} which represent the lower and upper truncation points respectively.
Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be carried out numerically.
}
\details{
The normal symmetric truncated distribution is a special case of the trucated normal distribution.
See \code{\link{Normal_trunc_ab}}.
}
\seealso{
\pkg{\link{ExtDist}} for other standard distributions.
}
\author{
Haizhen Wu and A. Jonathan R. Godfrey.
}
