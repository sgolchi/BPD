#### Preparing data for inference #####

## signal histograms from which signal functin is constructed
signal = read.table("data/HiggsSignalHistos.txt", sep="")

## simulated data (background + signal)
sim.data = read.table("data/data.dat", sep="")

## grid over mass range
G = signal[1:322,2]
n = length(G)

## signal function
sighat = function(G, m, s, c) {
  out = c * dnorm(G, m, s)
  return(out)
}

## estimating equation for signal function parameters
est.eq = function(par, signal, m) {
  s = par[1]
  c = par[2]
  shat = c * dnorm(G, m, s)
  out = sum((shat - signal) ^ 2)
  return(out)
}

## estimate signal function parameters
s = array(dim = c(3, 41))
c = array(dim = c(3, 41))
m = c(120, 125, 130)
for (j in 1:3) {
  mj = m[j]
  indexj = (41 * 322 * (j - 1)):(41 * 322 * j)
  for (i in 1:41) {
    indexi = (322 * (i - 1) + 1):(322 * i)
    opt = nlminb(c(2, .005), est.eq, lower = c(0, 0), signal = signal[indexj,][indexi,1], m = mj)
    par = opt$par
    s[j,i] = par[1]
    c[j,i] = par[2]
  }
}

## Extrapolate to find the signal functions for all the 322 masses (bin centres):
X = cbind(rep(1, 3), c(120, 125, 130), c(120 ^ 2, 125 ^ 2, 130 ^ 2))
betas = solve(X) %*% s
betac = solve(X) %*% c

Xnew = cbind(rep(1, 322), G, G ^ 2)
snew = Xnew %*% betas
cnew = Xnew %*% betac
snew[snew < 0] = 1e-6
cnew[cnew < 0] = 0
cnew[snew<=1e-6]<-0


break.vec = c(signal[1:322,2] - .25, signal[322,2])
bindata = hist(sim.data[,2], breaks = break.vec, plot = F)
Data = cbind(bindata$mids, bindata$counts)


Data[1:2,2] = Data[3,2] ## correcting edge effects
y = Data[,2]

H_data = list(grid = G, y = y, cnew = cnew, snew = snew)
save(H_data, file = 'data/Higgs_data')
