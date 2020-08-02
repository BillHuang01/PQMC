# load library
source("scripts/lib.R")
library(mvtnorm)
library(randtoolbox)

logmixture <- function(x){
  logf <- rep(0, 3)
  v <- (c(-5, 6, 3) + 20) / 40
  sigma <- 0.2
  for (i in 1:3){
    logf[i] <- sum(dnorm(x, mean = v[i], sd = sigma, log = T))
  }
  logf <- logaddexp(logf) -log(3)
  return (logf)
}

p <- 20
runs <- 100
steps <- 10
N <- 50
J <- 40
sigma <- 0.2

for (d in 2:p){
  print(d)
  set.seed(950922)
  expectation <- 0.5 + rep(4/3, d) / 40
  ini <- sobol(N, d)
  # store results
  mn.m.std <- rep(0,runs)
  mn.m.wts <- rep(0,runs)
  mn.z.std <- rep(0,runs)
  mn.z.wts <- rep(0,runs)
  mn.times <- rep(0,runs)
  ss.m.std <- rep(0,runs)
  ss.m.wts <- rep(0,runs)
  ss.z.std <- rep(0,runs)
  ss.z.wts <- rep(0,runs)
  ss.times <- rep(0,runs)
  sp.m.std <- rep(0,runs)
  sp.m.wts <- rep(0,runs)
  sp.m.las <- rep(0,runs)
  sp.z.std <- rep(0,runs)
  sp.z.wts <- rep(0,runs)
  sp.times <- rep(0,runs)
  for (k in 1:runs){
    # multinomial
    print("multinomial")
    start.time <- Sys.time()
    pmc.mn <- pmc(ini, logmixture, J, steps, sigma, resample = "Multinomial",
                  sigma.adapt = T, qmc = T, visualization = F)
    mn.times[k] <- as.numeric((Sys.time() - start.time), units = "secs")
    mn.m.std[k] <- log(mean((pmc.mn$m.std - expectation)^2))
    mn.m.wts[k] <- log(mean((pmc.mn$m.wts - expectation)^2))
    mn.z.std[k] <- log((pmc.mn$z.std - 1)^2)
    mn.z.wts[k] <- log((pmc.mn$z.wts - 1)^2)
    
    # systematic
    print("systematic")
    start.time <- Sys.time()
    pmc.ss <- pmc(ini, logmixture, J, steps, sigma, resample = "Systematic",
                  sigma.adapt = T, qmc = T,  visualization = F)
    ss.times[k] <- as.numeric((Sys.time() - start.time), units = "secs")
    ss.m.std[k] <- log(mean((pmc.ss$m.std - expectation)^2))
    ss.m.wts[k] <- log(mean((pmc.ss$m.wts - expectation)^2))
    ss.z.std[k] <- log((pmc.ss$z.std - 1)^2)
    ss.z.wts[k] <- log((pmc.ss$z.wts - 1)^2)
    
    # sp
    print("sp")
    start.time <- Sys.time()
    pmc.sp <- pmc(ini, logmixture, J, steps, sigma, resample = "SP",
                  sigma.adapt = T, qmc = T, visualization = F)
    sp.times[k] <- as.numeric((Sys.time() - start.time), units = "secs")
    sp.m.std[k] <- log(mean((pmc.sp$m.std - expectation)^2))
    sp.m.wts[k] <- log(mean((pmc.sp$m.wts - expectation)^2))
    sp.m.las[k] <- log(mean((pmc.sp$m.las - expectation)^2))
    sp.z.std[k] <- log((pmc.sp$z.std - 1)^2)
    sp.z.wts[k] <- log((pmc.sp$z.wts - 1)^2)
    
    # save log
    sink("results/pmc_hd/log.txt")
    cat(sprintf("N: %d\n", N))
    cat(sprintf("J: %d\n", J))
    cat(sprintf("P: %d\n", d))
    cat(sprintf("sigma: %.1f\n", sigma))
    cat(sprintf("run: %d\n", k))
    sink()
  }
  # save output
  logmse <- data.frame(
    mn.m.std,
    mn.m.wts,
    mn.z.std,
    mn.z.wts,
    mn.times,
    ss.m.std,
    ss.m.wts,
    ss.z.std,
    ss.z.wts,
    ss.times,
    sp.m.std,
    sp.m.wts,
    sp.m.las,
    sp.z.std,
    sp.z.wts,
    sp.times
  )
  file <- sprintf("results/pmc_hd/pmc_%dd_%d_%d_%d_%.1f_full_adapt.csv", d, N, J, steps, sigma)
  write.csv(logmse, file, row.names = F)
}
