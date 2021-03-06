% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calib_funs.R
\name{discovery_thresh}
\alias{discovery_thresh}
\title{Discovery threshold function}
\usage{
discovery_thresh(p0, G, snew, cnew, covmat, lambda_mean, alpha1, M)
}
\arguments{
\item{p0}{prior mass on \code{m = 0}}

\item{G}{grid over which the signal is evaluated}

\item{snew}{matrix of signal widths, same dimensions as \code{cnew}}

\item{cnew}{matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)}

\item{covmat}{covariance matrix for teh background prior}

\item{lambda_mean}{mean vector for the background prior}

\item{M}{number of iterations used to estimate the lower 5 percent quantile of the mass marginal posterior}

\item{lambda0}{initial backgorund value}

\item{y}{observations}
}
\value{
discovery threshold
}
\description{
Discovery threshold function
}
\examples{
load("data/SMC_full_poly")
L = dim(post2$thetat)[2]
dd_theta = density(post2$thetat[,L])
theta = dd_theta$x[which.max(dd_theta$y)]
dd_tau = density(post2$taut[,L])
tau = dd_tau$x[which.max(dd_tau$y)]
beta = apply(post2$beta[,L,], 2, mean)
xG = (G  - min(G))/(max(G) - min(G))
covmat = Gausscor(xG, theta, tau)
lambda_mean = log(poly(beta, xG)) - tau/2
cg = floor(seq(1, length(G), length = 100))
Qdis = discovery_thresh(p0 = .5, G, snew, cnew, covmat, lambda_mean, alpha1 = 3e-07, 1000)
}
