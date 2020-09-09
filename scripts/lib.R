# load library
library(pracma)
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

pmc <- function(logf, N, J, steps, ini, ini.logq = NULL,
                      sampling = "random", resampling = "sp", output = 'estimator',
                      sigma = NULL, sigma.adapt = T, visualization = F){
  # Population Monte Carlo with Prior Information
  # logf: log target density
  # ini: initialization points with log prior density,
  #      size of ini = N * J
  # N: number of proposal centers
  # J: number of samples per proposal
  # steps: number of PMC iterations
  # sample: random / qmc / sp / msp (mixture sp)
  # resample: multinomial / residual / systematic / stratified / sp
  # output: weighted samples / estimator
  # sigma: sigma for covariance, if NULL, then covariance adaptation
  # sigma.adapt: covariance adaptation ?
  # visualization: for two dimensional visualization
  
  # input size check
  if (nrow(ini) == N * J){
    ini.type <- 'sample'
    if (any(is.null(ini.logq))) stop('ini.logq cannot be null for sample ini!')
    if (length(ini.logq) != nrow(ini)) stop('ini.logq length does not match with ini #rows!')
  } else if (nrow(ini) == N){
    ini.type <- 'center'
    center <- ini
    if (is.null(sigma)) stop("if ini is proposal centers, then sigma is required!")
  } else {
    stop("number of ini points must be equal to N (for proposal centers) or N*J (for initial samples)!")
  }
  
  # visualization check
  p <- ncol(ini)
  if (visualization & p != 2) visualization <- F
  if (visualization){
    x1 <- x2 <- seq(0, 1, length.out = 101)
    x.grid <- expand.grid(x1, x2)
    z <- matrix(exp(apply(x.grid, 1, logf)), 101, 101)
  }
  
  # sample check
  if (sampling == "sp"){
    invisible(capture.output(noise.base <- sp(J, p, dist.str = rep("normal",p))$sp))
  }
  
  # information storage
  samp.all <- NULL
  samp.all.logwts <- NULL
  center.all <- NULL
  ess <- rep(NA, steps)

  if (ini.type == 'center') center.all <- rbind(center.all, center)
  
  # pmc procedures
  for (t in 1:steps){
    
    # sampling
    if (t == 1 & ini.type == 'sample'){
      samp <- ini # samples
      samp.logq <- ini.logq # log proposal density
    } else {
      # samples
      if (sampling == "random"){
        samp <- NULL
        for (i in 1:N){
          noise <- matrix(rnorm(J*p, mean = 0, sd = sigma), ncol = p)
          samp <- rbind(samp, (rep(1,J) %*% t(center[i,]) + noise))
        }
      } else if (sampling == 'qmc'){
        samp <- NULL
        for (i in 1:N){
          noise <- qnorm(sobol(J,p,scrambling=1,seed=sample(1e6,1)), mean = 0, sd = sigma)
          samp <- rbind(samp, (rep(1,J) %*% t(center[i,]) + noise))
        }
      } else if (sampling == "sp"){
        samp <- NULL
        for (i in 1:N){
          rotation <- randortho(p, type = "orthonormal")
          noise <- noise.base %*% t(rotation)
          noise <- noise[,sample(1:p,p)]
          noise <- noise * sigma
          samp <- rbind(samp, (rep(1,J) %*% t(center[i,]) + noise))
        }
      } else if (sampling == 'msp'){
        samp.qmc <- NULL
        for (i in 1:N){
          noise <- qnorm(sobol(J*10,p,scrambling=1,seed=sample(1e6,1)), mean = 0, sd = sigma)
          samp.qmc <- rbind(samp.qmc, (rep(1,J*10) %*% t(center[i,]) + noise))
        }
        invisible(capture.output(samp <- sp((N*J), p, dist.samp = samp.qmc)$sp))
        samp.nn.idx <- knnx.index(data = samp.qmc, query = samp, k = 1)
        samp <- samp.qmc[samp.nn.idx,]
      } else {
        stop (sprintf("invalid sampling method: %s", sampling))
      }
      # compute logq (deterministic mixture weight)
      samp.dm.logwts <- matrix(0, nrow = nrow(samp), ncol = N)
      for (i in 1:N){
        for (j in 1:nrow(samp)){
          dist <- samp[j,] - center[i,]
          samp.dm.logwts[j,i] <- sum(dnorm(dist, mean = 0, sd = sigma, log = T))
        }
      }
      samp.logq <- apply(samp.dm.logwts, 1, logaddexp) - log(N)
    }
    
    # weighting step
    samp.logf <- apply(samp, 1, logf)
    samp.logf[is.na(samp.logf)] <- -Inf
    samp.logwts <- samp.logf - samp.logq
    samp.wts <- exp(samp.logwts - logaddexp(samp.logwts))
    
    # store information
    samp.all <- rbind(samp.all, samp)
    samp.all.logwts <- c(samp.all.logwts, samp.logwts)
    ess[t] <- 1 / sum(samp.wts^2)
    # print(ess[t])
    
    # adaptation
    # covariance adaptation
    if (sigma.adapt & t < steps & !(t == 1 & ini.type == 'sample')){
      samp.dm.wts <- exp(samp.dm.logwts - (samp.logq + log(N)) %*% t(rep(1,N)))
      sigma <- 0
      for (i in 1:N){
        for (j in 1:nrow(samp)){
          sigma <- sigma + samp.wts[j] * samp.dm.wts[j,i] * sum((samp[j,] - center[i,])^2)
        }
      }
      sigma <- sqrt(sigma/p)
    }
      
    # center adaptation
    if (resampling == 'multinomial'){
      center <- samp[sample(1:length(samp.wts), N, replace = T, prob = samp.wts),]
    } else if (resampling == 'residual'){
      center <- samp[rs.sample(1:length(samp.wts), N, samp.wts),]
    } else if (resampling == 'systematic'){
      center <- samp[ss.sample(1:length(samp.wts), N, samp.wts),]
    } else if (resampling == 'stratified'){
      center <- samp[st.sample(1:length(samp.wts), N, samp.wts),]
    } else if (resampling == 'sp'){
      center <- samp[sp.sample(samp, N, samp.wts),]
    } else {
      stop (sprintf("invalid resampling method: %s", resampling))
    }
      
    center.all <- rbind(center.all, center)
      
    # sigma for initial step
    if (t == 1 & is.null(sigma)){
      center.dist <- as.matrix(dist(center))
      diag(center.dist) <- NA
      sigma <- max(apply(center.dist, 1, min, na.rm = T)) / sqrt(p)
      if (sigma == 0){
        center.dist <- rep(1,N*J) %*% t(center[1,]) - samp
        center.dist <- apply(center.dist, 1, function(x) sqrt(sum(x^2)))
        center.dist[center.dist == 0] <- NA
        sigma <- min(center.dist, na.rm = T) / sqrt(p)
      }
    }
    # print(sigma)
    
    # visulization
    if (visualization){
      contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 15, main = sprintf("t = %d", t))
      points(samp, pch = 18, cex = 0.5, col = "green")
      points(center, pch = 16, cex = 1, col = "red")
    }
    
  }
  
  if (output == 'estimator'){
    # standard pmc
    z.std <- mean(exp(samp.all.logwts))
    samp.std.wts <- exp(samp.all.logwts - logaddexp(samp.all.logwts))
    m.std <- c(t(samp.all) %*% samp.std.wts)
    # weighted pmc
    alpha <- ess / sum(ess)
    samp.wts.logwts <- samp.all.logwts + rep(log(alpha), each = N*J) + log(steps)
    z.wts <- mean(exp(samp.wts.logwts))
    # normalized as a whole
    samp.wts.wts <- exp(samp.wts.logwts - logaddexp(samp.wts.logwts))
    m.wts <- c(t(samp.all) %*% samp.wts.wts)
    # last adapted sample
    m.las <- c(apply(center, 2, mean))
    
    return (list(m.std = m.std,
                 m.wts = m.wts,
                 m.las = m.las,
                 z.std = z.std,
                 z.wts = z.wts))
  } else {
    return (list(samp.all = samp.all,
                 samp.all.logwts = samp.all.logwts,
                 center.all = center.all,
                 ess = ess))
  }
  
}
