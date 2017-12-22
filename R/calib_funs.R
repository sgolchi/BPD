

#' Signal function
#'
#' @param m the mass of Higgs
#' @param G grid over which the signal is evaluated
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @export
#'

sigfun = function(m, G, snew, cnew) {
  if (m == 0) sig = rep(0, length(G)) else {
    J = which.min(abs(G - m))
    SG = apply(cbind(cnew[J,], snew[J,]), 1, function(x) x[1] * dnorm(G, m, x[2]))
    sig = apply(SG, 1, sum)
  }
  return(sig)
}


bvec = function(lambda, sig, y) {
  term1 = (lambda - 1) * exp(lambda)
  term2 = (y * exp(lambda) / (exp(lambda) + sig)) * (1 - (sig * lambda / (exp(lambda) + sig)))
  return(term1 + term2)
}

cvec = function(lambda, sig, y) return(exp(lambda) * (1 - (y * sig) / (exp(lambda) + sig)^2))


lapapprox = function(m, lambda0, G, snew, cnew, y, covmat, lambda_mean) {
  sig = sigfun(m, G, snew, cnew)
  epsilon = 1e-6
  lambda = matrix(lambda0, 1, length(lambda0))
  t = 0
  repeat {
    t = t + 1
    b = bvec(lambda[t,], sig, y)
    c = cvec(lambda[t,], sig, y)
    mat = diag(1/c) + covmat
    matinv = solve(mat)
    A = covmat - covmat %*% matinv %*% covmat
    mean = ginv(covmat) %*% matrix(lambda_mean, length(lambda_mean), 1) + b
    lambda = rbind(lambda, t(A %*% mean))
    if (sum((lambda[t+1,] - lambda[t,])^2) <= epsilon) break
  }
  tl = nrow(lambda)
  lam_mean = lambda[tl, ] + diag(A)
  A = (A + t(A)) / 2
  lprm = dmnorm(lambda[tl,], lambda_mean, covmat, log = T) + log(.5/length(G))
  gam = exp(lam_mean) + sig
  llm = sum(dpois(y, gam, log = T))
  lpom = lprm + llm - dmnorm(lambda[tl,], lambda[tl,], A, log = T)
  return(lpom)
}


#' Mass marginal posterior density with the Laplace approximation
#'
#' @param j index
#' @param m the mass of Higgs
#' @param lambda0 initial backgorund value
#' @param G grid over which the signal is evaluated
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param y observations
#' @param covmat covariance matrix for teh background prior
#' @param lambda_mean mean vector for the background prior
#' @export
#'

simpost = function(j, m, G, snew, cnew, covmat, lambda_mean){
  lambda_sim = rmnorm(1, lambda_mean, covmat)
  sig = sigfun(m, G, snew, cnew)
  gam = exp(lambda_sim) + sig
  repeat {
    y<-c(); for (i in 1:length(G)) y[i] = rpois(1, gam[i])
    if (max(y)<5000) break
  }
  mass_rg = c(0, G)
  lpm = sapply(mass_rg, lapapprox, lambda0 = lambda_mean, G = G, snew = snew, cnew = cnew, y = y,
              covmat = covmat, lambda_mean= lambda_mean)
  lpm = lpm - max(lpm)
  pm = exp(lpm)
  pm = pm/sum(pm)
  J = which.min(abs(G - m))
  return(pm[J+1])
}

#' Mass marginal posterior density with the Laplace approximation with importance weights under "no Higgs" assumption
#'
#' @param j index
#' @param p0 prior mass on \code{m = 0}
#' @param lambda0 initial backgorund value
#' @param G grid over which the signal is evaluated
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param y observations
#' @param covmat covariance matrix for teh background prior
#' @param lambda_mean mean vector for the background prior
#' @export
#'

simpost_IS = function(j, p0, G, snew, cnew, covmat, lambda_mean){
  m = runif(1, 100, 180)
  lambda_sim = rmnorm(1, lambda_mean, covmat)
  sig = sigfun(m, G, snew, cnew)
  gam = exp(lambda_sim) + sig
  repeat {
    y<-c(); for (i in 1:length(G)) y[i] = rpois(1, gam[i])
    if (max(y)<5000) break
  }
  mass_rg = c(0, G)
  lpm = sapply(mass_rg, lapapprox, lambda0 = lambda_mean, G = G, snew = snew, cnew = cnew, y = y,
               covmat = covmat, lambda_mean= lambda_mean)
  lpm = lpm - max(lpm)
  pm = exp(lpm)
  pm = pm/sum(pm)
  W = ((1-p0) * pm[1])/(p0*(length(G)*sum(pm[2:(length(G)+1)])))
  J = which.min(abs(G - m))
  out = data.frame(pm = pm[1], W = W)
  return(out)
}

#' Exclusion threshold function
#'
#' @param m the mass of Higgs
#' @param lambda0 initial backgorund value
#' @param G grid over which the signal is evaluated
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param y observations
#' @param covmat covariance matrix for teh background prior
#' @param lambda_mean mean vector for the background prior
#' @param M number of iterations used to estimate the lower 5% quantile of the mass marginal posterior
#' @export
#' @return exclusion threshold for a given mass
#' @examples
#' load('data/SMC_full_poly')
#' L = dim(post2$thetat)[2]
#' dd_theta = density(post2$thetat[,L])
#' theta = dd_theta$x[which.max(dd_theta$y)]
#' dd_tau = density(post2$taut[,L])
#' tau = dd_tau$x[which.max(dd_tau$y)]
#' beta = apply(post2$beta[,L,], 2, mean)
#' xG = (G  - min(G))/(max(G) - min(G))
#' covmat = Gausscor(xG, theta, tau)
#' lambda_mean = log(poly(beta, xG)) - tau/2
#' cg = floor(seq(1, length(G), length = 100))
#' Qex = sapply(seq(100,180,1), exclusion_thresh,  G[cg], snew, cnew, covmat[cg,cg], lambda_mean[cg], M = 1000)

exclusion_thresh = function(m, G, snew, cnew, covmat, lambda_mean, alpha2, M) {
  pm = sapply(1:M, simpost, m = m, G = G, snew = snew, cnew = cnew, covmat = covmat,
             lambda_mean = lambda_mean)
  return(quantile(pm, p = alpha2))
}


calib_fun = function(q0, alpha, p, W){
  check = p < q0
  out = sum(p[check] * W[check]) - alpha
  return(out)
}

#' Discovery threshold function
#'
#' @param p0 prior mass on \code{m = 0}
#' @param lambda0 initial backgorund value
#' @param G grid over which the signal is evaluated
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param y observations
#' @param covmat covariance matrix for teh background prior
#' @param lambda_mean mean vector for the background prior
#' @param M number of iterations used to estimate the lower 5% quantile of the mass marginal posterior
#' @export
#' @return discovery threshold
#' @examples
#' load("data/SMC_full_poly")
#' L = dim(post2$thetat)[2]
#' dd_theta = density(post2$thetat[,L])
#' theta = dd_theta$x[which.max(dd_theta$y)]
#' dd_tau = density(post2$taut[,L])
#' tau = dd_tau$x[which.max(dd_tau$y)]
#' beta = apply(post2$beta[,L,], 2, mean)
#' xG = (G  - min(G))/(max(G) - min(G))
#' covmat = Gausscor(xG, theta, tau)
#' lambda_mean = log(poly(beta, xG)) - tau/2
#' cg = floor(seq(1, length(G), length = 100))
#' Qdis = discovery_thresh(p0 = .5, G, snew, cnew, covmat, lambda_mean, alpha1 = 3e-07, 1000)


discovery_thresh = function(p0, G, snew, cnew, covmat, lambda_mean, alpha1, M) {
  pw = t(sapply(1:M, simpost_IS, p0 = p0, G = G, snew = snew, cnew = cnew, covmat = covmat,
              lambda_mean = lambda_mean))
  q0 = uniroot(calib_fun, interval = c(min(unlist(pw[,1])), max(unlist(pw[,1]))), alpha = alpha1,
               p = unlist(pw[,1]), W = unlist(pw[,2]))$root
  return(q0)
}
