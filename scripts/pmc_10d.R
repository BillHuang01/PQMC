# load library
source("scripts/lib.R")
# source("lib.R")
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

# experiment setting
set.seed(950922)
p <- 10
expectation <- 0.5 + rep(4/3, p) / 40
Z <- 1
runs <- 100
steps <- 10
N <- 50 
J <- 40 
sampling <- "qmc" # random / qmc / sp / msp
resampling <- "sp" # multinomial / residual / systematic / stratified / sp

# initialization
ini <- sobol(N, p)
ini.logq <- NULL
ini.label <- "center"
#ini <- halton(N*J, p)
#ini.logq <- rep(0, N*J)
#ini.label <- "sample"

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
for (k in 1:runs){
  start.time <- Sys.time()
  pmc.output <- pmc(logmixture, N, J, steps, ini, ini.logq, 
                    sampling = sampling, resampling = resampling, 
                    sigma = sigma, sigma.adapt = sigma.adapt, visualization = F)
  end.time <- Sys.time()
  times[k] <- as.numeric((end.time - start.time), units = "secs")
  m.std[k] <- log(mean((pmc.output$m.std - expectation)^2))
  m.wts[k] <- log(mean((pmc.output$m.wts - expectation)^2))
  m.las[k] <- log(mean((pmc.output$m.las - expectation)^2))
  z.std[k] <- log((pmc.output$z.std - Z)^2)
  z.wts[k] <- log((pmc.output$z.wts - Z)^2)
  
  # save log
  log.file <- sprintf("results/pmc_10d/log_mixture_%d_%d_%d_%s_%s_%s_%s.txt",
                      N,J,steps,sampling,resampling,ini.label,sigma.label)
  sink(log.file)
  cat(sprintf("run: %d\n", k))
  cat(sprintf("m.std: %.3f\n", m.std[k]))
  cat(sprintf("m.wts: %.3f\n", m.wts[k]))
  cat(sprintf("time: %.3f\n", times[k]))
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
file <- sprintf("results/pmc_10d/pmc_mixture_%d_%d_%d_%s_%s_%s_%s.csv", 
                N, J, steps, sampling, resampling, ini.label, sigma.label)
write.csv(logmse, file, row.names = F)