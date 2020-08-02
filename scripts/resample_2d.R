source("scripts/lib.R")
library(mvtnorm)
library(randtoolbox)

# compare different resample techniques by visualization
logmixture <- function(x){
  #x <- (x - 0.5) * 40
  v <- matrix(NA, nrow = 5, ncol = 2)
  v[1,] <- c(-10,-10)
  v[2,] <- c(0,16)
  v[3,] <- c(13,8)
  v[4,] <- c(-9,7)
  v[5,] <- c(14,-14)
  v <- (v + 20) / 40
  sigma <- array(NA, dim = c(5,2,2))
  sigma[1,,] <- matrix(c(2,0.6,0.6,1), nrow = 2) / 40^2
  sigma[2,,] <- matrix(c(2,-0.4,-0.4,2), nrow = 2) / 40^2
  sigma[3,,] <- matrix(c(2,0.8,0.8,2), nrow = 2) / 40^2
  sigma[4,,] <- matrix(c(3,0,0,0.5), nrow = 2) / 40^2
  sigma[5,,] <- matrix(c(2,-0.1,-0.1,2), nrow = 2) / 40^2
  logf <- rep(0,5)
  for (i in 1:5) logf[i] <- dmvnorm(x, mean = v[i,], sigma = sigma[i,,], log = T)
  logf <- logaddexp(logf) - log(5)
  return (logf)
}

# contour information
x1 <- x2 <- seq(0, 1, length.out = 101)
x.grid <- expand.grid(x1, x2)
z <- matrix(exp(apply(x.grid, 1, logmixture)), 101, 101)

# importance samples
# importance support point
is.samp <- sobol(10000, 2)
is.logwts <- apply(is.samp, 1, logmixture)
is.wts <- exp(is.logwts - max(is.logwts))
is.wts <- is.wts / sum(is.wts)

# parameter
n <- 100
p <- 2
set.seed(950922)

# multinomial resampling
mn.samp <- is.samp[sample(1:10000, n, replace = T, prob = is.wts),]
contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 10)
points(mn.samp, pch = 16, cex = 1, col = "red")

# residual resampling
rs.samp <- is.samp[rs.sample(1:10000, n, is.wts),]
contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 10)
points(rs.samp, pch = 16, cex = 1, col = "red")

# stratified resampling
st.samp <- is.samp[st.sample(1:10000, n, is.wts),]
contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 10)
points(st.samp, pch = 16, cex = 1, col = "red")

# systematic resampling
ss.samp <- is.samp[ss.sample(1:10000, n, is.wts),]
contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 10)
points(ss.samp, pch = 16, cex = 1, col = "red")

# support points resampling
sp.samp <- is.samp[sp.sample(is.samp, n, is.wts),]
contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 10)
points(sp.samp, pch = 16, cex = 1, col = "red")
