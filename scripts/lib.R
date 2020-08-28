# load library
library(support)
library(FNN)

# log to handle numerical underflow
logaddexp <- function(logv){
  logv.max <- max(logv)
  logv.sum <- log(sum(exp(logv - logv.max))) + logv.max
  return(logv.sum)
}

# resample function
rs.sample <- function(x, n, prob){
  # residual sampling
  idx <- c()
  N <- length(x)
  ept <- n * prob
  cnt <- floor(ept)
  wts <- ept - cnt
  if (sum(wts) > 0){
    wts <- wts / sum(wts)
    idx <- sample(1:N, (n-sum(cnt)), replace = T, prob = wts)
  }
  for (i in 1:N) idx <- c(idx, rep(i,cnt[i]))
  return (x[idx])
}

ss.sample <- function(x, n, prob){
  # systematic sampling
  idx <- c()
  N <- length(x)
  ridx <- sample(1:N, N, replace = F)
  x <- x[ridx]
  prob <- prob[ridx]
  u <- runif(1) / n
  l <- 0
  j <- 0
  while (u < 1){
    if (l > u){
      u <- u + 1 / n
      idx <- c(idx, j)
    } else {
      j <- j + 1
      l <- l + prob[j]
    }
  }
  return(x[idx])
}

st.sample <- function(x, n, prob){
  # stratified sampling
  idx <- c()
  N <- length(x)
  ridx <- sample(1:N, N, replace = F)
  x <- x[ridx]
  prob <- prob[ridx]
  prob.cumsum <- cumsum(prob)
  endpoints <- seq(0, 1, length.out = n + 1)
  for (i in 1:n){
    lb <- sum(prob.cumsum <= endpoints[i]) + 1
    ub <- min(sum(prob.cumsum <= endpoints[(i+1)]) + 1, N)
    if (lb == ub){
      idx <- c(idx, lb)
    } else {
      prob.loc <- prob[c(lb:ub)]
      prob.loc[1] <- prob.cumsum[lb] - endpoints[i]
      prob.loc[length(prob.loc)] <- endpoints[(i+1)] - prob.cumsum[(ub-1)]
      prob.loc <- prob.loc / sum(prob.loc)
      idx <- c(idx, sample(lb:ub, 1, replace = T, prob = prob.loc))
    }
  }
  return (x[idx])
}

sp.sample <- function(x, n, prob, tol = 1e-6, iter.max = 10){
  # support points
  if (is.null(dim(x))) stop("x must be a matrix!")
  N <- nrow(x)
  x.dist <- as.matrix(dist(x))
  x.measure <- c(x.dist %*% prob)
  # initial sample
  idx <- c(which.min(x.measure))
  if (prob[idx[1]] < 1/N) idx[1] <- which.max(prob)
  for (i in 2:n){
    measure <- x.measure - apply(matrix(x.dist[,idx],ncol=(i-1)), 1, sum) / i
    idx <- c(idx, which.min(measure))
  }
  edist <- 2 * sum(x.measure[idx]) / n - sum(x.dist[idx,idx]) / n^2
  # iterative update for improvement
  iter <- 0
  while (TRUE){
    iter <- iter + 1
    for (i in 1:n){
      measure <- x.measure - apply(x.dist[,idx[-i]], 1, sum) / n
      idx[i] <- which.min(measure)
    }
    edist.new <- 2 * sum(x.measure[idx]) / n - sum(x.dist[idx,idx]) / n^2
    if ((edist-edist.new) < tol || iter > iter.max){
      break
    } else {
      edist <- edist.new
    }
  }
  return (idx)
}

pmc <- function(ini, logf, J, steps, sigma, 
                resample = 'SP', sample = "random", output = 'estimator',
                sigma.adapt = F, visualization = F){
  # Population Monte Carlo Global Resampling
  # ini: initialization points
  # logf: log target density
  # J: number of samples per proposal
  # sigma: initial proposal standard deviation
  # resample: Multinomial/Residual/Systematic/SP
  # output: esimation / raw
  # sigma.adapt: adaption for sigma (T/F)
  # qmc: qmc proposal points (T/F)
  # visualization: T/F
  
  if (visualization){
    x1 <- x2 <- seq(0, 1, length.out = 101)
    x.grid <- expand.grid(x1, x2)
    z <- matrix(exp(apply(x.grid, 1, logf)), 101, 101)
  }
  
  n <- nrow(ini)
  p <- ncol(ini)
  # store samples
  ini.all <- ini
  samp.all <- NULL
  samp.all.logwts <- NULL
  samp.no <- rep(0, steps)
  ess <- rep(0, steps)
  
  for (t in 1:steps){
    
    # sample from proposals
    if (sample == 'random'){
      samp <- NULL
      for (i in 1:n){
        noise <- matrix(rnorm(J*p, mean = 0, sd = sigma), ncol = p)
        samp <- rbind(samp, (rep(1,J) %*% t(ini[i,]) + noise))
      }
    } else if (sample == 'qmc'){
      samp <- NULL
      for (i in 1:n){
        noise <- qnorm(sobol(J,p,scrambling=1,seed=sample(1e6,1)), mean = 0, sd = sigma)
        samp <- rbind(samp, (rep(1,J) %*% t(ini[i,]) + noise))
      }
    } else if (sample == 'sp'){
      samp.qmc <- NULL
      for (i in 1:n){
        noise <- qnorm(sobol(J*25,p,scrambling=1,seed=sample(1e6,1)), mean = 0, sd = sigma)
        samp.qmc <- rbind(samp.qmc, (rep(1,J*25) %*% t(ini[i,]) + noise))
      }
      invisible(capture.output(samp <- sp((n*J), p, dist.samp = samp.qmc)$sp))
      samp.nn.idx <- knnx.index(data = samp.qmc, query = samp, k = 1)
      samp <- samp.qmc[samp.nn.idx,]
    }
    
    # compute deterministic mixture weight
    samp.dm.logwts <- matrix(0, nrow = nrow(samp), ncol = n)
    for (i in 1:n){
      for (j in 1:nrow(samp)){
        dist <- samp[j,] - ini[i,]
        samp.dm.logwts[j,i] <- sum(dnorm(dist, mean = 0, sd = sigma, log = T))
      }
    }
    samp.dm.logsumwts <- apply(samp.dm.logwts, 1, logaddexp)
    samp.logwts <- apply(samp, 1, logf) - (samp.dm.logsumwts - log(n))
    samp.wts <- exp(samp.logwts - logaddexp(samp.logwts))
    
    samp.all <- rbind(samp.all, samp)
    samp.all.logwts <- c(samp.all.logwts, samp.logwts)
    ess[t] <- 1 / sum(samp.wts^2)
    
    if (sigma.adapt & t < steps){
      samp.dm.wts <- exp(samp.dm.logwts - samp.dm.logsumwts %*% t(rep(1,n)))
      sigma <- 0
      for (i in 1:n){
        for (j in 1:nrow(samp)){
          sigma <- sigma + samp.wts[j] * samp.dm.wts[j,i] * sum((samp[j,] - ini[i,])^2)
        }
      }
      sigma <- sqrt(sigma/p)
    }
    
    # resample: multinomial/residual/systematic/sp
    if (resample == 'multinomial'){
      ini <- samp[sample(1:length(samp.wts), n, replace = T, prob = samp.wts),]
    } else if (resample == 'residual'){
      ini <- samp[rs.sample(1:length(samp.wts), n, samp.wts),]
    } else if (resample == 'systematic'){
      ini <- samp[ss.sample(1:length(samp.wts), n, samp.wts),]
    } else if (resample == 'stratified'){
      ini <- samp[st.sample(1:length(samp.wts), n, samp.wts),]
    } else if (resample == 'sp'){
      ini <- samp[sp.sample(samp, n, samp.wts),]
    } else {
      stop("no such resampling method!")
    }
    
    ini.all <- rbind(ini.all, ini)
    
    if (visualization){
      contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 15, main = sprintf("t = %d", t))
      points(samp, pch = 18, cex = 0.5, col = "green")
      points(ini, pch = 16, cex = 1, col = "red")
    }
    
  }
  
  if (output == 'estimator'){
    # standard pmc
    z.std <- mean(exp(samp.all.logwts))
    samp.std.wts <- exp(samp.all.logwts - logaddexp(samp.all.logwts))
    m.std <- c(t(samp.all) %*% samp.std.wts)
    # weighted pmc
    alpha <- ess / sum(ess)
    samp.wts.logwts <- samp.all.logwts + rep(log(alpha), each = n*J) + log(steps)
    z.wts <- mean(exp(samp.wts.logwts))
    samp.wts.wts <- exp(samp.wts.logwts - logaddexp(samp.wts.logwts))
    m.wts <- c(t(samp.all) %*% samp.wts.wts)
    # last adapted sample
    m.las <- c(apply(ini, 2, mean))
    
    return (list(m.std = m.std,
                 m.wts = m.wts,
                 m.las = m.las,
                 z.std = z.std,
                 z.wts = z.wts))
  } else {
    return (list(samp.all = samp.all,
                 samp.all.logwts = samp.all.logwts,
                 ini.all = ini.all,
                 ess = ess))
  }
  
}

energy_dist <- function(X, Y, X.wts = NULL, Y.wts = NULL){
  n.X <- nrow(X)
  n.Y <- nrow(Y)
  if (is.null(X.wts)) X.wts <- rep(1/n.X, n.X)
  if (is.null(Y.wts)) Y.wts <- rep(1/n.Y, n.Y)
  e.dist <- 0
  for (i in 1:n.X){
    for (j in 1:n.Y){
      e.dist <- e.dist + 2 * X.wts[i] * Y.wts[j] * sqrt(sum((X[i,]-Y[j,])^2))
    }
  }
  for (i in 1:n.X){
    for (j in 1:n.X){
      e.dist <- e.dist - X.wts[i] * X.wts[j] * sqrt(sum((X[i,]-X[j,])^2))
    }
  }
  #for (i in 1:n.Y){
  #   for (j in 1:n.Y){
  #      e.dist <- e.dist - Y.wts[i] * Y.wts[j] * sqrt(sum((Y[i,]-Y[j,])^2))
  #   }
  #}
  return (e.dist)
}