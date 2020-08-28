# load library
source("scripts/lib.R")
# source("lib.R")
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
runs <- 2
steps <- 8
p <- 2
N <- 25 # 25, 50
J <- 10 # 10, 5
sample <- "sp" # random / qmc / sp
resample <- "sp" # multinomial / systematic / sp
ini <- sobol(N, p)
ini.label <- "full"
# ini <- 0.4 + 0.2 * sobol(N, p)
# ini.label <- "sub"

stds <- c(0.1,0.2,0.5)
if (sample == "sp" | resample == "sp"){
  adaptation <- c(TRUE)
  adaptation.label <- c("adapt")
} else {
  adaptation <- c(TRUE, FALSE)
  adaptation.label <- c("adapt","regular")
}

for (i in 1:length(adaptation)){
  for (j in 1:length(stds)){
    sigma <- stds[j]
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
      pmc.output <- pmc(ini, logmixture, J, steps, sigma, 
                        sample = sample, resample = resample, 
                        sigma.adapt = adaptation[i], visualization = F)
      end.time <- Sys.time()
      times[k] <- as.numeric((end.time - start.time), units = "secs")
      m.std[k] <- log(mean((pmc.output$m.std - expectation)^2))
      m.wts[k] <- log(mean((pmc.output$m.wts - expectation)^2))
      m.las[k] <- log(mean((pmc.output$m.las - expectation)^2))
      z.std[k] <- log((pmc.output$z.std - 1)^2)
      z.wts[k] <- log((pmc.output$z.wts - 1)^2)
      
      # save log
      log.file <- sprintf("results/pmc_2d/log_%d_%s_%s_%s.txt",
                          N,sample,resample,ini.label)
      sink(log.file)
      cat(sprintf("N: %d\n", N))
      cat(sprintf("J: %d\n", J))
      cat(sprintf("sigma: %.1f\n", sigma))
      cat(sprintf("adaptation: %s\n", adaptation.label[i]))
      cat(sprintf("run: %d\n", k))
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
    file <- sprintf("results/pmc_2d/pmc_2d_%d_%d_%d_%.1f_%s_%s_%s_%s.csv", 
                    N, J, steps, sigma, sample, resample, ini.label, adaptation.label[i])
    write.csv(logmse, file, row.names = F)
  }
}