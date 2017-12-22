


#' Bernstein Polynomial
#'
#' @param beta polynomila coefficients
#' @param x the point at which the polynomial is evaluated
#' @param n order of the polynomial
#' @export
#'

poly = function(beta, x, n = 4) {
  b = array(dim = c(5, length(x)))
  for (i in 0:n) {
    b[i + 1,] = dbinom(i, n, x)
  }
  f = matrix(beta, 1, 5) %*% b
  return(f)
}


#' Gaussian covariance function
#'
#' @param x vector of points for which the covariance matrix is generated
#' @param theta correlation parameter
#' @param tau variance parameter
#' @export
#' @return A covariance matrix of dimension ,equal to the length of \code{x}


Gausscor = function(x, theta, tau){
  n  = length(x)
  Sigma_0  = array(dim = c(n, n));
  for (i in 1:n) {
    for (j in 1:n) {
      Sigma_0[i,j] = tau * exp(- theta * (x[i] - x[j]) ^ 2 )
    }
  }
  return(Sigma_0 + diag(rep(.0001, n)))
}

#' Tempered log likelihood function
#'
#' @param par vector of parameters including background, mass and signal strength
#' @param y vector of observations
#' @param G grid over mass range
#' @param temp temperature parameter
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @export
#' @return the tempered log likelihood for a set of parameters

temploglik = function(par, y = y, G = G, temp, cnew, snew) {
  n = length(y)
  lambda = par[1:n]
  mass = par[n + 1]
  mu = par[n + 2]
  if (mass==0) gam = lambda else {
    J = which.min(abs(G - mass))
    SG = apply(cbind(cnew[J,], snew[J,]), 1, function(x) x[1] * dnorm(G, mass, x[2]))
    SG = apply(SG, 1, sum)
    gam = lambda + mu * SG
  }
  loglik = sum(dpois(y, gam, log = TRUE))
  out = temp * loglik
  return(out)
}


#' Log posterior density function
#'
#' @param lambda vector of background values, the same size of grid \code{G}
#' @param mass mass value
#' @param mu_s signal streangth parameter
#' @param theta correlation parameter of the background logGP prior
#' @param tau variance parameter of the background logGP prior
#' @param beta polynomial coefficients of the background mean of the background logGP prior
#' @param y vector of observations
#' @param G grid over mass range
#' @param temp temperature parameter
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param p0 prior probability of $mass(Higgs)=0$
#' @param mu_var prior variance of signal strength parameter
#' @param back_mean parametric form used to specify the mean of background logGP prior
#' @export
#' @return the value of log posterior density for a set of parameters

logpost = function(lambda, mass, mu_s, theta, tau, beta, y = y, G = G, temp, cnew, snew,
                   p0 = .5, mu_var, back_mean = 'B polynomial'){
  x = (G - min(G)) / (max(G) - min(G))
  Rz = Gausscor(x, theta, tau) #+ diag(rep(1e-3, 322))
  if (back_mean == 'exponential decay') {
    muG = exp_decay(x, beta)
    lp3 = dnorm(beta[1], 1000, 100, log = T) + sum(dnorm(beta[2], 3, 1), log = T)
  }
  if (back_mean == 'B polynomial') {
    muG = poly(beta, x, n = 4)
    lp3 = dnorm(beta[1], 1000, 100, log = T) + sum(dnorm(beta[2:5], 200, 50), log = T)
  }
  if (back_mean == 'constant') {
    muG = rep(beta, length(y))
    lp3 = dnorm(beta, mean(y), 100, log = T)
  }
  muz = c(log(muG) - tau / 2)
  lp0 = (mass == 0) * log(p0) + (mass > 0) * log((1 - p0) / (180 - 100))
  lp1 = -log(tau) + dgamma(theta, 1, 1, log = T) #+ log(dinvgamma(tau, 1, 1))
  lp2 = dnorm(log(mu_s), 0, mu_var, log = T)
  logprior = lp0 + lp1 + lp2 + lp3 + dmnorm(log(lambda), muz, Rz, log = T)
  - sum(log(lambda))
  if (mass == 0) gam = lambda else {
    J = which.min(abs(G - mass))
    SG = apply(cbind(cnew[J,], snew[J,]), 1, function(x) x[1] * dnorm(G, mass, x[2]))
    SG = apply(SG, 1, sum)
    gam = lambda + mu_s * SG
  }
  loglik = sum(dpois(y, gam, log = T))
  out = logprior + temp * loglik
  return(out)
}


#' Adaptive tempering helper function
#'
#' @param temp current temperature
#' @param temp0 previous temperature
#' @param N Monte Carlo sample size
#' @param mt current sample of mass, vector of size \code{N}
#' @param mut current sample of signal strength parameter, vector of size \code{N}
#' @param lambdat current sample of baclground paths, array of dimensions \code{N*length(G)}
#' @param G grid over mass range
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param y vector of observations
#' @export
#' @return current temperature to attain an effective sample size of \code{N/2}

adapt_seq = function(temp, temp0, N, mt, mut, lambdat, G, cnew, snew, y) {
  logden = apply(cbind(lambdat, mt, mut), 1, temploglik, y = y,
                 G = G, temp = temp0, cnew, snew)
  lognum = apply(cbind(lambdat, mt, mut), 1, temploglik, y = y,
                 G = G, temp = temp, cnew, snew)
  w = exp(lognum - logden)
  w = w / sum(w, na.rm = T)
  ESS = ifelse(sum(is.na(w)) == N, 0, 1 / sum(w ^ 2, na.rm = T))
  out = ESS - (N / 2)
  if (temp == temp0) out = N/2
  return(out)
}


#' SMC sampler
#'
#' @param N Monte Carlo sample size
#' @param G grid over mass range
#' @param y vector of observations
#' @param cnew matrix of signal normalizing constants - dimensions are (size of grid \code{G})by(number of categories/channels,41)
#' @param snew matrix of signal widths, same dimensions as \code{cnew}
#' @param p0 prior probability of mass(Higgs)=0
#' @param mu_var prior variance of signal strength parameter
#' @param back_mean parametric form used to specify the mean of background logGP prior
#' @export
#' @return A list of sequence of posterior samples for model parameters,
#' sequence of acceptance rates for all parameters (for monitoring),
#' sequence of effective sample sizes (ESS) and sequence of temperatures
#' @examples
#' load('data/Higgs_data')
#' G = H_data$grid
#' y = H_data$y
#' cnew = H_data$cnew
#' snew = H_data$snew
#' post = SMC(N = 1e03, G = G, y = y, cnew = cnew, snew = snew, p0 = .5, mu_var = .5, back_mean = 'B polynomial')
#' save(post, file = 'data/SMC_full_poly')

SMC = function(N = 1e03, G, y, cnew, snew, p0 = .5, mu_var = .5,
               back_mean = 'B polynomial') {
  tempseq = 0
  tempseq_T = 1
  n = length(y)
  x = (G - min(G)) / (max(G) - min(G))
  lpdent = array(NA, dim = c(N, 1))
  lambdat = array(NA, dim = c(N, 1, n))
  ESS = c()
  pi0 = c()
  atheta = 0
  atau = 0
  amu = 0
  am = 0
  alambda = 0
  qa = .01
  q1 = .1
  t = 1
  W = array(rep(1 / N, N), dim = c(N, 1))
  if (back_mean == 'exponential decay') {
    beta = array(NA, dim = c(N, 1, 2))
    beta[,t,1] = rnorm(N, 1000, 100)
    beta[,t,2] = rnorm(N, 3, 1)
    abeta = array(rep(0, 2), dim = c(2, 1))
  }
  if (back_mean == 'B polynomial') {
    beta = array(NA, dim = c(N, 1, 5))
    beta[,t,1] = rnorm(N, 1000, 100)
    beta[,t,2:5] = matrix(rnorm(4 * N, 200, 50), N, 4)
    abeta = array(rep(0, 5), dim = c(5, 1))
  }
  if (back_mean == 'constant') {
    beta = array(NA, dim = c(N, 1, 1))
    beta[,t,] = rnorm(N, mean(y), 100)
    abeta = array(rep(0, 1), dim = c(1, 1))
  }
  mt = array(NA, dim = c(N, 1))
  for (i in 1:N) {
    u = runif(1)
    if (u <= .5) mt[i,t] = 0
    else mt[i,t] = runif(1, 100, 180)
  }
  thetat = array(rgamma(N, 1, 1), dim = c(N, 1))
  taut = array((rnorm(N, 0, 10))^2, dim = c(N, 1)) #rinvgamma(N, 1, 1)
  mu_st = array(exp(rnorm(N, 0, mu_var)), dim = c(N, 1))
  bad = which(abs(mu_st[,t] - 1) < 1e-6)
  mu_st[bad,1] = mu_st[bad,1] + 1e-6
  for (j in 1:N) {
    Rz = Gausscor(x, thetat[j,t], taut[j,t])# + diag(rep(1e-3, 322))
    if (back_mean == 'exponential decay') {
      muG = exp_decay(x, beta[j,t,])
    }
    if (back_mean == 'B polynomial') {
      muG = poly(beta[j,t,], x, n = 4)
    }
    if (back_mean == 'constant') {
      muG = rep(beta[j,t], length(y))
    }
    muz = c(log(muG) - taut[j,t] / 2)
    lambdat[j,t,] =  exp(rmnorm(1, muz, Rz))
    lpdent[j,t] = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], thetat[j,t], taut[j,t],
                          beta[j,t,], y, G, temp = tempseq[1], cnew, snew, p0, mu_var, back_mean)
  }
  qq = sd(mt[, 1])
  pi0[1] = length(mt[mt[,1] == 0, 1]) / N
  den = density(mt[mt[,1] > 0, 1], n = 1000, from = 100, to = 180)
  kernx = den$x
  kerny = den$y
  repeat {
    t = t + 1
    qa = c(qa, .01 / t)
    atheta = c(atheta, 0)
    atau = c(atau, 0)
    amu = c(amu, 0)
    am = c(am, 0)
    alambda = c(alambda, 0)
    abeta = abind(abeta, rep(0, length(abeta[,t-1])), along = 2)
    lpdent = abind(lpdent, lpdent[,t - 1], along = 2)
    thetat = abind(thetat, thetat[,t - 1], along = 2)
    taut = abind(taut, taut[,t - 1], along = 2)
    mu_st = abind(mu_st, mu_st[,t - 1], along = 2)
    beta = abind(beta, beta[,t-1,], along = 2)
    lambdat = abind(lambdat, lambdat[,t - 1,], along = 2)
    mt = abind(mt, mt[,t - 1], along = 2)
    W = abind(W, W[,t - 1], along = 2)
    tempseq[t] = ifelse(adapt_seq(temp = tempseq_T, temp0 = tempseq[t-1], N = N, mt = mt[,t],
                                  mut = mu_st[,t], lambdat[,t,],G = G, cnew, snew, y = y) > 0,
                       tempseq_T, uniroot(adapt_seq, interval = c(tempseq_T, tempseq[t-1]),
                                         temp0 = tempseq[t-1], N = N, mt = mt[,t], mut = mu_st[,t],
                                         lambdat[,t,], G = G, cnew, snew, y = y)$root)
    logden = apply(cbind(lambdat[,t,], mt[,t], mu_st[,t]), 1, temploglik, y = y,
                   G = G, temp = tempseq[t - 1], cnew, snew)
    lognum = apply(cbind(lambdat[,t,], mt[,t], mu_st[,t]), 1, temploglik, y = y,
                   G = G, temp = tempseq[t], cnew, snew)
    w = exp(lognum - logden)
    w = w / sum(w, na.rm = T)
    ESS[t] = ifelse(sum(is.na(w)) == N, 0, 1 / sum(w ^ 2, na.rm = T))
    index = sample(1:N, prob = w, replace = TRUE)
    lpdent[,t] = lpdent[index,t]
    thetat[,t] = thetat[index,t]
    taut[,t] = taut[index,t]
    mu_st[,t] = mu_st[index,t]
    mt[,t] = mt[index,t]
    lambdat[,t,] = lambdat[index,t,]
    beta[,t,] = beta[index,t,]
    #W[,t] = rep(1 / N, N)
    qq = sd(mt[,t])
    qmu = sd(mu_st[,t])
    qbeta = apply(beta[,t,], 2, sd)
    for (j in 1:N) {
      ## sampling theta
      repeat {
        newtheta  = rchisq(1, thetat[j,t])
        if (newtheta>.005) break
      }
      qnum = dchisq(thetat[j,t], newtheta)
      lpnum = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], newtheta, taut[j,t], beta[j,t,],
                      y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
      qden = dchisq(newtheta, thetat[j,t])
      lpdent[j,t] = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], thetat[i,t], taut[j,t],
                            beta[j,t,], y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
      p = min(1, exp(lpnum + qnum - lpdent[j,t] - qden))
      if (runif(1) <= p) {
        thetat[j,t] = newtheta
        atheta[t] = atheta[t] + 1
        lpdent[j,t] = lpnum
      }
      ## sampling tau
      repeat{
        newtau = rchisq(1, taut[j,t])
        if (newtau>.001) break
      }
      qnum = dchisq(taut[j,t], newtau, log = T)
      qden = dchisq(newtau, taut[j,t], log = T)
      lpnum = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], thetat[j,t], newtau,
                      beta[j,t,], y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
      p = min(1, exp(lpnum + qnum - lpdent[j,t] - qden))
      if (runif(1) <= p) {
        taut[j,t] = newtau
        atau[t] = atau[t] + 1
        lpdent[j,t] = lpnum
      }

      ## sampling beta
      if (back_mean == 'exponential decay') {
        for (i in 1:2) {
          newbeta = beta[j,t,]
          betai = rnorm(1, beta[j,t,i], qbeta[i] / 10)
          newbeta[i] = betai
          lpnum = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], thetat[j,t], taut[j,t], newbeta, y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
          p = min(1, exp(lpnum - lpdent[j,t]))
          if (runif(1) <= p) {
            beta[j,t,i] = betai
            abeta[i,t] = abeta[i,t] + 1
            lpdent[j,t] = lpnum
          }
        }
      }
      if (back_mean == 'B polynomial') {
        for (i in 1:5) {
          newbeta = beta[j,t,]
          betai = rnorm(1, beta[j,t,i], qbeta[i] / 10)
          newbeta[i] = betai
          lpnum = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], thetat[j,t], taut[j,t], newbeta, y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
          p = min(1, exp(lpnum - lpdent[j,t]))
          if (runif(1) <= p) {
            beta[j,t,i] = betai
            abeta[i,t] = abeta[i,t] + 1
            lpdent[j,t] = lpnum
          }
        }
      }
      if (back_mean == 'constant') {
        newbeta = rnorm(1, beta[j,t,], qbeta / 10)
        lpnum = logpost(lambdat[j,t,], mt[j,t], mu_st[j,t], thetat[j,t], taut[j,t], newbeta, y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
        p = min(1, exp(lpnum - lpdent[j,t]))
        if (runif(1) <= p) {
          beta[j,t] = newbeta
          abeta[t] = abeta[t] + 1
          lpdent[j,t] = lpnum
        }
      }
      ## sampling lambda
      Rz = Gausscor(x, thetat[j,t], taut[j,t]) #+ diag(rep(1e-3, 322))
      prop.cov = qa[t] * Rz
      new0 = rmnorm(1, log(lambdat[j,t,]), prop.cov)
      qnum = - sum(log(lambdat[j,t,]))
      lpnum = logpost(exp(new0), mt[j,t], mu_st[j,t], thetat[j,t], taut[j,t],
                      beta[j,t,], y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
      qden = - sum(new0)
      p = min(1, exp(lpnum - lpdent[j,t] + qnum - qden))
      if (runif(1) <= p) {
        lambdat[j,t,] = exp(new0)
        alambda[t] = alambda[t] + 1
        lpdent[j,t] = lpnum
      }
      ## sampling mu
      lognewmu = rnorm(1, log(mu_st[j,t]), qmu)
      qnum = - log(mu_st[j,t])
      lpnum = logpost(lambdat[j,t,], mt[j,t], exp(lognewmu), thetat[j,t], taut[j,t],
                      beta[j,t,], y = y, G = G, temp = tempseq[t], cnew, snew, p0, mu_var, back_mean)
      qden = - lognewmu
      p = min(1, exp(lpnum - lpdent[j,t] + qnum - qden))
      if (runif(1) <= p) {
        mu_st[j,t] = exp(lognewmu)
        amu[t] = amu[t] + 1
        lpdent[j,t] = lpnum
      }
      ## sampling mass
      v = runif(1)
      if (v <= pi0[t - 1]) mnew = 0 else mnew = sample(kernx, 1, prob = kerny)
      lpnum = logpost(lambdat[j,t,], mnew, mu_st[j,t], thetat[j,t], taut[j,t],
                      beta[j,t,], y = y, G = G, temp = tempseq[t], cnew, snew, p0 = p0, mu_var = mu_var, back_mean)
      if (mt[j,t] == 0) qnum = log(pi0[t - 1]) else qnum = log(kerny[which(abs(kernx - mt[j,t]) < .05)])
      if (mnew == 0) qden = log(pi0[t - 1]) else qden = log(kerny[which(abs(kernx - mnew) < .05)])
      p = min(1, exp(lpnum + qnum - lpdent[j,t] - qden))
      if (runif(1) <= p) {
        mt[j,t] = mnew
        am[t] = am[t] + 1
        lpdent[j,t] = lpnum
      }
    }
    pi0[t] = length(mt[mt[,t] == 0,t]) / N
    den = density(mt[mt[,t] > 0, t], n = 1000, from = 100, to = 180)
    kernx = den$x
    kerny = den$y
    if (tempseq[t] >= tempseq_T) break
  }
  out = list(thetat = thetat, taut = taut, mt = mt, mu_st = mu_st, beta = beta,
             lambdat = lambdat, atheta = atheta, atau = atau, am = am, amu = amu,
             alambda = alambda, abeta = abeta, ESS = ESS, tempseq = tempseq)
  return(out)
}
