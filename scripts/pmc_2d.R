# load library
#source("scripts/lib.R")
source("lib.R")
library(mvtnorm)
library(randtoolbox)

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

# experiment setting
set.seed(950922)
expectation <- 0.5 + c(1.6,1.4) / 40
runs <- 100
steps <- 10
p <- 2
N <- 50 # 25, 50, 100
J <- 20 # 40, 20, 10
stds <- c(0.1,0.2,0.5)
ini <- sobol(N, p)
ini.label <- "full"
# ini <- 0.4 + 0.2 * sobol(N, p)
# ini.label <- "sub"

adaptation <- c(TRUE, FALSE)
adaptation.label <- c("adapt","regular")
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
                    sigma.adapt = adaptation[i], qmc = T, visualization = T)
      mn.times[k] <- Sys.time() - start.time
      mn.m.std[k] <- log(mean((pmc.mn$m.std - expectation)^2))
      mn.m.wts[k] <- log(mean((pmc.mn$m.wts - expectation)^2))
      mn.z.std[k] <- log((pmc.mn$z.std - 1)^2)
      mn.z.wts[k] <- log((pmc.mn$z.wts - 1)^2)
      
      # systematic
      print("systematic")
      start.time <- Sys.time()
      pmc.ss <- pmc(ini, logmixture, J, steps, sigma, resample = "Systematic",
                    sigma.adapt = adaptation[i], qmc = T,  visualization = T)
      ss.times[k] <- Sys.time() - start.time
      ss.m.std[k] <- log(mean((pmc.ss$m.std - expectation)^2))
      ss.m.wts[k] <- log(mean((pmc.ss$m.wts - expectation)^2))
      ss.z.std[k] <- log((pmc.ss$z.std - 1)^2)
      ss.z.wts[k] <- log((pmc.ss$z.wts - 1)^2)
      
      # sp
      print("sp")
      start.time <- Sys.time()
      pmc.sp <- pmc(ini, logmixture, J, steps, sigma, resample = "SP",
                    sigma.adapt = adaptation[i], qmc = T, visualization = T)
      sp.times[k] <- Sys.time() - start.time
      sp.m.std[k] <- log(mean((pmc.sp$m.std - expectation)^2))
      sp.m.wts[k] <- log(mean((pmc.sp$m.wts - expectation)^2))
      sp.m.las[k] <- log(mean((pmc.sp$m.las - expectation)^2))
      sp.z.std[k] <- log((pmc.sp$z.std - 1)^2)
      sp.z.wts[k] <- log((pmc.sp$z.wts - 1)^2)
      
      # save log
      sink("results/pmc_2d/log.txt")
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
    file <- sprintf("results/pmc_2d/pmc_2d_%d_%d_%d_%.1f_%s_%s.csv", 
                    N, J, steps, sigma, ini.label, adaptation.label[i])
    write.csv(logmse, file, row.names = F)
  }
}