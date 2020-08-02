source("scripts/lib.R")
library(randtoolbox)

# same sigma for all dimensions
set.seed(950922)
N <- 1000
n <- 100
p <- 20
sigma <- sqrt(2)
runs <- 100

ess <- rep(0, p)
is.logmse <- rep(0, p)
sp.ep.logmse <- rep(0, p)
sp.is.logmse <- rep(0, p)
sp.times <- rep(0, p)
for (d in 2:p){
  print(d)
  samp <- matrix(qnorm(sobol(N,d),0,sigma), ncol = d)
  samp.logf <- apply(samp, 1, function(x) sum(dnorm(x,0,1,log=T)))
  samp.logq <- apply(samp, 1, function(x) sum(dnorm(x,0,sigma,log=T)))
  samp.logwts <- samp.logf - samp.logq
  samp.wts <- exp(samp.logwts - max(samp.logwts))
  samp.wts <- samp.wts / sum(samp.wts)
  ess[d] <- 1 / sum(samp.wts^2)
  is.est <- c(t(samp) %*% samp.wts)
  is.logmse[d] <- log(mean(is.est^2))
  # resample
  mn.ep.logmse <- rep(0, runs)
  mn.is.logmse <- rep(0, runs)
  mn.times <- rep(0, runs)
  rs.ep.logmse <- rep(0, runs)
  rs.is.logmse <- rep(0, runs)
  rs.times <- rep(0, runs)
  st.ep.logmse <- rep(0, runs)
  st.is.logmse <- rep(0, runs)
  st.times <- rep(0, runs)
  ss.ep.logmse <- rep(0, runs)
  ss.is.logmse <- rep(0, runs)
  ss.times <- rep(0, runs)
  for (j in 1:runs){
    # multinomial
    start.time <- Sys.time()
    samp.mn <- samp[sample(1:N, n, replace = T, prob = samp.wts),]
    mn.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    mn.est <- apply(samp.mn, 2, mean)
    mn.ep.logmse[j] <- log(mean(mn.est^2))
    mn.is.logmse[j] <- log(mean((mn.est - is.est)^2))
    # residual
    start.time <- Sys.time()
    samp.rs <- samp[rs.sample(1:N, n, samp.wts),]
    rs.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    rs.est <- apply(samp.rs, 2, mean)
    rs.ep.logmse[j] <- log(mean(rs.est^2))
    rs.is.logmse[j] <- log(mean((rs.est - is.est)^2))
    # systematic
    start.time <- Sys.time()
    samp.st <- samp[st.sample(1:N, n, samp.wts),]
    st.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    st.est <- apply(samp.st, 2, mean)
    st.ep.logmse[j] <- log(mean(st.est^2))
    st.is.logmse[j] <- log(mean((st.est - is.est)^2))
    # stratified
    start.time <- Sys.time()
    samp.ss <- samp[ss.sample(1:N, n, samp.wts),]
    ss.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    ss.est <- apply(samp.ss, 2, mean)
    ss.ep.logmse[j] <- log(mean(ss.est^2))
    ss.is.logmse[j] <- log(mean((ss.est - is.est)^2))
  }
  # store output
  logmse <- data.frame(
    mn.ep = mn.ep.logmse,
    mn.is = mn.is.logmse,
    mn.times = mn.times,
    rs.ep = rs.ep.logmse,
    rs.is = rs.is.logmse,
    rs.times = mn.times,
    st.ep = st.ep.logmse,
    st.is = st.is.logmse,
    st.times = st.times,
    ss.ep = ss.ep.logmse,
    ss.is = ss.is.logmse,
    ss.times = ss.times
  )
  file <- sprintf("results/resample/resample_ess_%dd.csv", d)
  write.csv(logmse, file, row.names = F)
  # sp resampling
  start.time <- Sys.time()
  samp.sp <- samp[sp.sample(samp, n, samp.wts),]
  sp.times[d] <- as.numeric((Sys.time() - start.time), units = "secs")
  sp.est <- apply(samp.sp, 2, mean)
  sp.ep.logmse[d] <- log(mean(sp.est^2))
  sp.is.logmse[d] <- log(mean((sp.est - is.est)^2))
}
logmse <- data.frame(
  dimension = c(1:p),
  ess = ess,
  is.logmse = is.logmse,
  sp.ep.logmse = sp.ep.logmse,
  sp.is.logmse = sp.is.logmse,
  sp.times = sp.times
)
file <- sprintf("results/resample/resample_ess_summary.csv")
write.csv(logmse, file, row.names = F)

# different sigma for different dimensions
set.seed(950922)
N <- 1000
n <- 100
p <- 20
runs <- 100

sigmas <- rep(0, p)
ess <- rep(0, p)
is.logmse <- rep(0, p)
sp.ep.logmse <- rep(0, p)
sp.is.logmse <- rep(0, p)
sp.times <- rep(0, p)
for (d in 2:p){
  print(d)
  sigma <- exp(2*log(3)/d^(0.8))
  sigmas[d] <- sigma
  samp <- matrix(qnorm(sobol(N,d),0,sigma), ncol = d)
  samp.logf <- apply(samp, 1, function(x) sum(dnorm(x,0,1,log=T)))
  samp.logq <- apply(samp, 1, function(x) sum(dnorm(x,0,sigma,log=T)))
  samp.logwts <- samp.logf - samp.logq
  samp.wts <- exp(samp.logwts - max(samp.logwts))
  samp.wts <- samp.wts / sum(samp.wts)
  ess[d] <- 1 / sum(samp.wts^2)
  is.est <- c(t(samp) %*% samp.wts)
  is.logmse[d] <- log(mean(is.est^2))
  # resample
  mn.ep.logmse <- rep(0, runs)
  mn.is.logmse <- rep(0, runs)
  mn.times <- rep(0, runs)
  rs.ep.logmse <- rep(0, runs)
  rs.is.logmse <- rep(0, runs)
  rs.times <- rep(0, runs)
  st.ep.logmse <- rep(0, runs)
  st.is.logmse <- rep(0, runs)
  st.times <- rep(0, runs)
  ss.ep.logmse <- rep(0, runs)
  ss.is.logmse <- rep(0, runs)
  ss.times <- rep(0, runs)
  for (j in 1:runs){
    # multinomial
    start.time <- Sys.time()
    samp.mn <- samp[sample(1:N, n, replace = T, prob = samp.wts),]
    mn.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    mn.est <- apply(samp.mn, 2, mean)
    mn.ep.logmse[j] <- log(mean(mn.est^2))
    mn.is.logmse[j] <- log(mean((mn.est - is.est)^2))
    # residual
    start.time <- Sys.time()
    samp.rs <- samp[rs.sample(1:N, n, samp.wts),]
    rs.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    rs.est <- apply(samp.rs, 2, mean)
    rs.ep.logmse[j] <- log(mean(rs.est^2))
    rs.is.logmse[j] <- log(mean((rs.est - is.est)^2))
    # systematic
    start.time <- Sys.time()
    samp.st <- samp[st.sample(1:N, n, samp.wts),]
    st.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    st.est <- apply(samp.st, 2, mean)
    st.ep.logmse[j] <- log(mean(st.est^2))
    st.is.logmse[j] <- log(mean((st.est - is.est)^2))
    # stratified
    start.time <- Sys.time()
    samp.ss <- samp[ss.sample(1:N, n, samp.wts),]
    ss.times[j] <- as.numeric((Sys.time() - start.time), units = "secs")
    ss.est <- apply(samp.ss, 2, mean)
    ss.ep.logmse[j] <- log(mean(ss.est^2))
    ss.is.logmse[j] <- log(mean((ss.est - is.est)^2))
  }
  # store output
  logmse <- data.frame(
    mn.ep = mn.ep.logmse,
    mn.is = mn.is.logmse,
    mn.times = mn.times,
    rs.ep = rs.ep.logmse,
    rs.is = rs.is.logmse,
    rs.times = mn.times,
    st.ep = st.ep.logmse,
    st.is = st.is.logmse,
    st.times = st.times,
    ss.ep = ss.ep.logmse,
    ss.is = ss.is.logmse,
    ss.times = ss.times
  )
  file <- sprintf("results/resample/resample_dimen_%dd.csv", d)
  write.csv(logmse, file, row.names = F)
  # sp resampling
  start.time <- Sys.time()
  samp.sp <- samp[sp.sample(samp, n, samp.wts),]
  sp.times[d] <- as.numeric((Sys.time() - start.time), units = "secs")
  sp.est <- apply(samp.sp, 2, mean)
  sp.ep.logmse[d] <- log(mean(sp.est^2))
  sp.is.logmse[d] <- log(mean((sp.est - is.est)^2))
}
logmse <- data.frame(
  dimension = c(1:p),
  sigma = sigmas,
  ess = ess,
  is.logmse = is.logmse,
  sp.ep.logmse = sp.ep.logmse,
  sp.is.logmse = sp.is.logmse,
  sp.times = sp.times
)
file <- sprintf("results/resample/resample_dimen_summary.csv")
write.csv(logmse, file, row.names = F)
