source("scripts/hilbert.R")

# log to handle numerical underflow
logaddexp <- function(logv){
  logv.max <- max(logv)
  logv.sum <- log(sum(exp(logv - logv.max))) + logv.max
  return(logv.sum)
}

# hilbert curve ordering
hilbert.curve.order <- function(x,order=8){
  dim <- ncol(x)
  for (i in 1:dim){
    xi.min <- min(x[,i])
    xi.max <- max(x[,i])
    xi.range <- xi.max - xi.min
    xi.min <- xi.min - 1/(2^(order-1)) * xi.range
    xi.max <- xi.max + 1/(2^(order-1)) * xi.range
    x[,i] <- as.integer((x[,i]-xi.min)/(xi.max-xi.min)*2^(order))
  }
  x.hdist <- apply(x, 1, hilbertcurve_distance_from_coordinates, dim=dim, order=order)
  x.order.idx <- order(x.hdist)
  return (x.order.idx)
}

# resample function
ms.sample <- function(x, n, prob){
  # multinomial resampling
  N <- nrow(x)
  idx <- sample(1:N, n, replace = T, prob = prob)
  return (idx)
}

rs.sample <- function(x, n, prob){
  # residual resampling
  idx <- c()
  N <- nrow(x)
  ept <- n * prob
  cnt <- floor(ept)
  wts <- ept - cnt
  if (sum(wts) > 0){
    wts <- wts / sum(wts)
    idx <- sample(1:N, (n-sum(cnt)), replace = T, prob = wts)
  }
  for (i in 1:N) idx <- c(idx, rep(i,cnt[i]))
  return (idx)
}

ss.sample <- function(x, n, prob){
  # systematic resampling
  N <- nrow(x)
  # order the samples by hilbert curve
  x.order.idx <- hilbert.curve.order(x)
  # reorder the probability vector
  prob <- prob[x.order.idx]
  # sampling
  idx <- c()
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
  return(x.order.idx[idx])
}

st.sample <- function(x, n, prob){
  # stratified resampling
  N <- nrow(x)
  # order the samples by hilbert curve
  x.order.idx <- hilbert.curve.order(x)
  # reorder the probability vector
  prob <- prob[x.order.idx]
  # sampling
  idx <- c()
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
  return (x.order.idx[idx])
}

sp.sample <- function(x, n, prob, tol = 1e-6, iter.max = 10){
  # support points resampling
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

pmc <- function(logf, K, J, steps, ini,
                sampling = "random", resampling = "sp", output = 'estimator',
                sigma = NULL, sigma.adapt = T, visualization = F){
  # Population Monte Carlo with Prior Information
  # logf: log target density
  # K: number of proposal centers
  # J: number of samples per proposal
  # ini: initialized centers
  # steps: number of PMC iterations
  # sample: random / qmc
  # resample: multinomial / residual / systematic / stratified / sp
  # output: weighted samples / estimator
  # sigma: sigma for covariance, if NULL, then covariance adaptation
  # sigma.adapt: covariance adaptation?
  # visualization: for two dimensional visualization
  
  # initialization check
  # center
  if (nrow(ini) == K){
    center <- ini
  } else {
    stop("number of initial centers must be equal to K!")
  }
  # sigma
  if (is.null(sigma)){
    if (sigma.adapt){
      sigma <- min(dist(center))
    } else {
      stop("sigma cannot be null if no adaptation!")
    }
  }
  
  # visualization check
  p <- ncol(ini)
  if (visualization & p != 2) visualization <- F
  if (visualization){
    x1 <- x2 <- seq(0, 1, length.out = 101)
    x.grid <- expand.grid(x1, x2)
    z <- matrix(exp(apply(x.grid, 1, logf)), 101, 101)
    # plot initial centers
    contour.default(x = x1, y = x2, z = z, drawlabels = F, nlevels = 15, main = sprintf("t = %d", 0))
    points(center, pch = 16, cex = 1, col = "red")
  }
  
  # information storage
  samp.all <- NULL
  samp.all.logwts <- NULL
  center.all <- NULL
  ess <- rep(NA, steps)
  center.all <- rbind(center.all, center)
  
  # pmc procedures
  for (t in 1:steps){
    # sampling
    if (sampling == "random"){
      samp <- NULL
      for (i in 1:K){
        noise <- matrix(rnorm(J*p, mean = 0, sd = sigma), ncol = p)
        samp <- rbind(samp, (rep(1,J) %*% t(center[i,]) + noise))
      }
    } else if (sampling == 'qmc'){
      samp <- NULL
      for (i in 1:K){
        noise <- qnorm(sobol(J,p,scrambling=1,seed=sample(1e6,1)), mean = 0, sd = sigma)
        samp <- rbind(samp, (rep(1,J) %*% t(center[i,]) + noise))
      }
    } else {
      stop (sprintf("invalid sampling method: %s", sampling))
    }
    # compute logq (deterministic mixture weight)
    samp.dm.logwts <- matrix(0, nrow = nrow(samp), ncol = K)
    for (i in 1:K){
      for (j in 1:nrow(samp)){
        dist <- samp[j,] - center[i,]
        samp.dm.logwts[j,i] <- sum(dnorm(dist, mean = 0, sd = sigma, log = T))
      }
    }
    samp.logq <- apply(samp.dm.logwts, 1, logaddexp) - log(K)
    
    # weighting step
    samp.logf <- apply(samp, 1, logf)
    samp.logf[is.na(samp.logf)] <- -Inf
    samp.logwts <- samp.logf - samp.logq
    samp.wts <- exp(samp.logwts - logaddexp(samp.logwts))
    
    # store information
    samp.all <- rbind(samp.all, samp)
    samp.all.logwts <- c(samp.all.logwts, samp.logwts)
    ess[t] <- 1 / sum(samp.wts^2)
    
    # adaptation
    # covariance adaptation
    if (sigma.adapt & t < steps){
      samp.dm.wts <- exp(samp.dm.logwts - (samp.logq + log(K)) %*% t(rep(1,K)))
      sigma <- 0
      for (i in 1:K){
        for (j in 1:nrow(samp)){
          sigma <- sigma + samp.wts[j] * samp.dm.wts[j,i] * sum((samp[j,] - center[i,])^2)
        }
      }
      sigma <- sqrt(sigma/p)
    }
      
    # center adaptation
    if (resampling == 'multinomial'){
      center <- samp[ms.sample(samp, K, samp.wts),]
    } else if (resampling == 'residual'){
      center <- samp[rs.sample(samp, K, samp.wts),]
    } else if (resampling == 'systematic'){
      center <- samp[ss.sample(samp, K, samp.wts),]
    } else if (resampling == 'stratified'){
      center <- samp[st.sample(samp, K, samp.wts),]
    } else if (resampling == 'sp'){
      center <- samp[sp.sample(samp, K, samp.wts),]
    } else {
      stop (sprintf("invalid resampling method: %s", resampling))
    }
      
    center.all <- rbind(center.all, center)
    
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
    samp.wts.logwts <- samp.all.logwts + rep(log(alpha), each = K*J) + log(steps)
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
