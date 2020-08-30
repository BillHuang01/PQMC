# load library
source("scripts/lib.R")
source("lib.R")
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

logmixture <- function(x){
  p <- length(x)
  v <- matrix(NA, nrow = 3, ncol = p)
  v[1,] <- rep(-5,p)
  v[2,] <- rep(3,p)
  v[3,] <- rep(8,p)
  v <- (20 + v) / 40
  sigma <- array(NA, dim = c(3,p,p))
  sigma.base <- diag(p)
  sigma[1,,] <- 0.15^2 * diag(p)
  sigma[2,,] <- 0.2^2 * (-0.5)^(abs(row(sigma.base) - col(sigma.base)))
  sigma[3,,] <- 0.2^2 * 0.5^(abs(row(sigma.base) - col(sigma.base)))
  logf <- rep(0,3)
  for (i in 1:3) logf[i] <- dmvnorm(x, mean = v[i,], sigma = sigma[i,,], log = T)
  logf <- logaddexp(logf) -log(3)
  return (logf)
}

# experiment setting
set.seed(950922)
p <- 10
expectation <- 0.5 + rep(4/3, p) / 40
runs <- 100
steps <- 8
N <- 50
J <- 40
sample <- "sp"
resample <- "sp"
ini <- sobol(N, p)
ini.label <- "full"

stds <- c(0.1,0.2,0.5)
if (sample == "sp" | resample == "sp"){
  adaptation <- c(TRUE)
  adaptation.label <- c("adapt")
} else {
  adaptation <- c(TRUE, FALSE)
  adaptation.label <- c("adapt","regular")
}

sigma <- 0.5
pmc.output <- pmc(ini, logmixture, J, steps, sigma, 
                  sample = sample, resample = resample,
                  sigma.adapt = T, visualization = F)

log(mean((pmc.output$m.std - expectation)^2))
log(mean((pmc.output$m.wts - expectation)^2))
log(mean((pmc.output$m.las - expectation)^2))
log((pmc.output$z.std - 1)^2)
log((pmc.output$z.wts - 1)^2)


for (i in 1:length(adaptation)){
  for (j in 1:length(stds)){
    sigma <- stds[j]
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
    # experiment
    for (k in 1:runs){
      # multinomial
      print("multinomial")
      start.time <- Sys.time()
      pmc.mn <- pmc(ini, logmixture, J, steps, sigma, resample = "Multinomial",
                    sigma.adapt = adaptation[i], qmc = T, visualization = F)
      mn.times[k] <- Sys.time() - start.time
      mn.m.std[k] <- log(mean((pmc.mn$m.std - expectation)^2))
      mn.m.wts[k] <- log(mean((pmc.mn$m.wts - expectation)^2))
      mn.z.std[k] <- log((pmc.mn$z.std - 1)^2)
      mn.z.wts[k] <- log((pmc.mn$z.wts - 1)^2)
      
      # systematic
      print("systematic")
      start.time <- Sys.time()
      pmc.ss <- pmc(ini, logmixture, J, steps, sigma, resample = "Systematic",
                    sigma.adapt = adaptation[i], qmc = T,  visualization = F)
      ss.times[k] <- Sys.time() - start.time
      ss.m.std[k] <- log(mean((pmc.ss$m.std - expectation)^2))
      ss.m.wts[k] <- log(mean((pmc.ss$m.wts - expectation)^2))
      ss.z.std[k] <- log((pmc.ss$z.std - 1)^2)
      ss.z.wts[k] <- log((pmc.ss$z.wts - 1)^2)
      
      # sp
      print("sp")
      start.time <- Sys.time()
      pmc.sp <- pmc(ini, logmixture, J, steps, sigma, resample = "SP",
                    sigma.adapt = adaptation[i], qmc = T, visualization = F)
      sp.times[k] <- Sys.time() - start.time
      sp.m.std[k] <- log(mean((pmc.sp$m.std - expectation)^2))
      sp.m.wts[k] <- log(mean((pmc.sp$m.wts - expectation)^2))
      sp.m.las[k] <- log(mean((pmc.sp$m.las - expectation)^2))
      sp.z.std[k] <- log((pmc.sp$z.std - 1)^2)
      sp.z.wts[k] <- log((pmc.sp$z.wts - 1)^2)
      
      # save log
      sink("results/pmc_10d/log.txt")
      cat(sprintf("N: %d\n", N))
      cat(sprintf("J: %d\n", J))
      cat(sprintf("sigma: %.1f\n", sigma))
      cat(sprintf("adaptation: %s\n", adaptation.label[i]))
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
    file <- sprintf("results/pmc_10d/pmc_10d_%d_%d_%d_%.1f_%s_%s.csv", 
                    N, J, steps, sigma, ini.label, adaptation.label[i])
    write.csv(logmse, file, row.names = F)
  }
}