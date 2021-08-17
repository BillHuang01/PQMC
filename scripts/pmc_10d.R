# read in arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("two arguments has to be provided!")

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
  logf <- logaddexp(logf) - log(3)
  return (logf)
}

# target distribution
logf <- logmixture
expectation <- 0.5 + rep(4/3, 10) / 40
Z <- 1
logf.label <- "mixture"

# experiment setting
set.seed(20210522)
runs <- 2
steps <- 10 # number of steps
p <- 10 # dimensions
K <- 50 # number of centers
J <- 40 # number of samples per center
sampling <- args[1] # random / qmc
resampling <- args[2] # multinomial / systematic / stratified / sp

# initialization
ini <- sobol(K, p)
# variance
sigma <- 0.2
sigma.adapt <- T
sigma.label <- "0.2_adapt"

# store results
m.std <- rep(0, runs)
m.wts <- rep(0, runs)
m.las <- rep(0, runs)
z.std <- rep(0, runs)
z.wts <- rep(0, runs)
times <- rep(0, runs)

# experiment
for (l in 1:runs){
  start.time <- Sys.time()
  pmc.output <- pmc(logmixture, K, J, steps, ini, 
                    sampling = sampling, resampling = resampling, 
                    sigma = sigma, sigma.adapt = sigma.adapt, visualization = F)
  end.time <- Sys.time()
  times[l] <- as.numeric((end.time - start.time), units = "secs")
  m.std[l] <- log(mean((pmc.output$m.std - expectation)^2))
  m.wts[l] <- log(mean((pmc.output$m.wts - expectation)^2))
  m.las[l] <- log(mean((pmc.output$m.las - expectation)^2))
  z.std[l] <- log((pmc.output$z.std - Z)^2)
  z.wts[l] <- log((pmc.output$z.wts - Z)^2)
  
  # save log
  log.file <- sprintf("results/pmc_10d/log_%s_%d_%d_%d_%s_%s_%s.txt",
                      logf.label,K,J,steps,sampling,resampling,sigma.label)
  sink(log.file)
  cat(sprintf("run: %d\n", l))
  cat(sprintf("m.std: %.3f\n", m.std[l]))
  cat(sprintf("m.wts: %.3f\n", m.wts[l]))
  cat(sprintf("time: %.3f\n", times[l]))
  sink()
}
# save output
logmse <- data.frame(
  m.std,
  m.wts,
  m.las,
  z.std,
  z.wts,
  times
)
file <- sprintf("results/pmc_10d/pmc_%s_%d_%d_%d_%s_%s_%s.csv", 
                logf.label, K, J, steps, sampling, resampling, sigma.label)
write.csv(logmse, file, row.names = F)