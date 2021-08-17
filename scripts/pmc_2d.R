# read in arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("two arguments has to be provided!")

# load library
source("scripts/lib.R")
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

logbanana <- function(x)
{
  if (any(x < 0)||any(x > 1)) return (log(0))
  lower1 <- -40
  upper1 <- 40
  lower2 <- -25
  upper2 <- 10
  theta1 <- lower1+(upper1-lower1)*x[1]
  theta2 <- lower2+(upper2-lower2)*x[2]+.03*theta1^2-3
  val <- -.5*(theta1^2/100+theta2^2)
  return(val)
}

# compute expectation by grid for banana shape distribution
# grid-based result
# library(cubature)
# banana <- function(x){
#   return (exp(logbanana(x)))
# }
# x1 <- x2 <- seq(0, 1, length = 1000)
# denom <- adaptIntegrate(banana, c(0,0),c(1,1))$int # 0.0223886
# f1 <- function(theta1) apply(cbind(theta1, theta2), 1, banana)
# f1.val <- rep(NA, 1000)
# for (i in 1:1000){
#   theta2 <- x2[i]
#   f1.val[i] <- integrate(f1, 0, 1)$val
# }
# f1.val <- f1.val / denom
# f1.val <- f1.val / sum(f1.val) * 1000
# mean(f1.val * x1) # 0.7162834
# f2 <- function(theta2) apply(cbind(theta1, theta2), 1, banana)
# f2.val <- rep(NA, 1000)
# for (i in 1:1000){
#   theta1 <- x1[i]
#   f2.val[i] <- integrate(f2, 0, 1)$val
# }
# f2.val <- f2.val / denom
# f2.val <- f2.val / sum(f2.val) * 1000
# mean(f2.val * x2) # 0.5

# target distribution
logf <- logmixture
expectation <- 0.5 + c(1.6,1.4) / 40
Z <- 1
logf.label <- "mixture"
# logf <- logbanana
# expectation <- c(0.5,0.7162834)
# z <- 0.0223886
# logf.label <- "banana"

# experiment setting
set.seed(20210522)
runs <- 100
steps <- 10 # number of steps
p <- 2 # dimensions
K <- 25 # number of centers
J <- 40 # number of samples per center
sampling <- args[1] # random / qmc
resampling <- args[2] # multinomial / residual / systematic / stratified / sp

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
  pmc.output <- pmc(logf, K, J, steps, ini, 
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
  log.file <- sprintf("results/pmc_2d/log_%s_%d_%d_%d_%s_%s_%s.txt",
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
file <- sprintf("results/pmc_2d/pmc_%s_%d_%d_%d_%s_%s_%s.csv", 
                logf.label, K, J, steps, sampling, resampling, sigma.label)
write.csv(logmse, file, row.names = F)