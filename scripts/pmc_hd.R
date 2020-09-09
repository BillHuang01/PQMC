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
  logf <- logaddexp(logf) -log(3)
  return (logf)
}

# experimental setting
p <- 15
expectation <- 0.5 + rep(4/3, p) / 40
runs <- 100
steps <- 10
N <- 50 
J <- 40
sampling <- "qmc" # random / qmc / sp / msp
resampling <- "sp" # multinomial / residual / systematic / stratified / sp

# initialization
ini.label <- "center"
# ini.label <- "sample"

# variance
sigma <- NULL
sigma.adapt <- T
sigma.label <- "adapt"

dir.create(sprintf("results/pmc_hd/mixture_%s_%s_%s_%s/", 
                   sampling, resampling, ini.label, sigma.label))

for (d in 2:p){
  set.seed(950922)
  expectation <- 0.5 + rep(4/3, d) / 40
  ini <- sobol(N, p)
  ini.logq <- NULL
  #ini <- halton(N*J, d)
  #ini.logq <- rep(0, N*J)
  
  # store results
  m.std <- rep(0, runs)
  m.wts1 <- rep(0, runs)
  m.wts2 <- rep(0, runs)
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
    m.wts1[k] <- log(mean((pmc.output$m.wts1 - expectation)^2))
    m.wts2[k] <- log(mean((pmc.output$m.wts2 - expectation)^2))
    m.las[k] <- log(mean((pmc.output$m.las - expectation)^2))
    z.std[k] <- log((pmc.output$z.std - 1)^2)
    z.wts[k] <- log((pmc.output$z.wts - 1)^2)
    
    # save log
    log.file <- sprintf("results/pmc_hd/log_mixture_%d_%d_%d_%s_%s_%s_%s.txt",
                        N,J,steps,sampling,resampling,ini.label,sigma.label)
    sink(log.file)
    cat(sprintf("dimension:%d\n", d))
    cat(sprintf("run: %d\n", k))
    cat(sprintf("m.std: %.3f\n", m.std[k]))
    cat(sprintf("m.wts1: %.3f\n", m.wts1[k]))
    cat(sprintf("m.wts2: %.3f\n", m.wts2[k]))
    cat(sprintf("time: %.3f\n", times[k]))
    sink()
  }
  
  # save output
  logmse <- data.frame(
    m.std,
    m.wts1,
    m.wts2,
    m.las,
    z.std,
    z.wts,
    times
  )
  file <- sprintf("results/pmc_hd/mixture_%s_%s_%s_%s/pmc_%dd_%d_%d_%d.txt", 
                  sampling, resampling, ini.label, sigma.label, d, N, J, steps)
  write.csv(logmse, file, row.names = F)
  
}