source("scripts/lib.R")
library(GPfit)
library(randtoolbox)

# log density function
# load data
fem <- read.csv("data/FEM.csv", header = T)
exp <- read.csv("data/EXP.csv", header = T)
# truncate x values to be no more than 11
fem <- fem[fem$x<11,]
exp <- exp[exp$x<11,]
# scale x to be between [0,1]
x.scale <- max(c(fem$x, exp$x))
exp.x <- exp$x / x.scale
exp.y <- exp$z

# load GP model
gp.model <- readRDS("data/drilling_gp.rds")

logprior <- function(x, a, b, lambda_a, lambda_b){
  val1 <- exp(lambda_a*(x-a)) * (x<a)
  val2 <- (x>=a) * (x<=b)
  val3 <- exp(-lambda_b*(x-b))  * (x>b)
  return(log(val1 + val2 + val3))
}

logf <- function(x){
  eta <- x[1]
  gamma1 <- x[2]
  gamma2 <- x[3]
  logeta <- logprior(eta, 0.5, 1, 10, 10)
  loggamma1 <- logprior(gamma1, 0.5, 1, 10, 100)
  loggamma2 <- logprior(gamma2, 0.75, 1.25, 10, 10)
  pred <- exp(apply(cbind(eta, gamma1*exp.x^gamma2), 1, 
                    function(x) predict.GP(gp.model, xnew=matrix(x,nrow=1))$Y_hat))
  logsse <- -length(exp.x)/2*log(sum((exp.y-pred)^2))
  return((logsse+logeta+loggamma1+loggamma2))
}


# experiment setting
set.seed(950922)
steps <- 7
p <- 3
N <- 13
J <- 7
sampling <- "qmc"
resampling <- "sp"

# initialization
ini <- as.matrix(read.table(sprintf("data/drilling_ini_%d.txt", N)))
ini.logq <- NULL
ini.label <- "center"
#ini <- halton(N*J, p)
#ini[,1] <- 0.5 + 0.5 * ini[,1]
#ini[,2] <- 0.5 + 0.5 * ini[,2]
#ini[,3] <- 0.75 + 0.5 * ini[,3]
#ini.logq <- rep(0, N*J)
#for (i in 1:(N*J)){
#  ini.logq[i] <- logprior(ini[i,1], 0.5, 1, 10, 10) + 
#    logprior(ini[i,2], 0.5, 1, 10, 100) + logprior(ini[i,3], 0.75, 1.25, 10, 10)
#}
#ini.label <- "sample"

# variance
sigma <- NULL
sigma.adapt <- T
sigma.label <- "adapt"

pmc.output <- pmc(logf, N, J, steps, ini, ini.logq, 
                  sampling = sampling, resampling = resampling, output = 'raw',
                  sigma = sigma, sigma.adapt = sigma.adapt, visualization = F)
output.label <- sprintf('%s_%s_%s_%s', sampling, resampling, ini.label, sigma.label)
write.table(pmc.output$samp.all, 
            sprintf("results/drilling/pmc_%s_samp.txt", output.label),
            row.names = F, col.names = F, sep = "\t")
write.table(matrix(pmc.output$samp.all.logwts, ncol = 1), 
            sprintf("results/drilling/pmc_%s_samp_logwts.txt", output.label),
            row.names = F, col.names = F)
write.table(pmc.output$center.all, 
            sprintf("results/drilling/pmc_%s_center.txt", output.label),
            row.names = F, col.names = F, sep = "\t")
write.table(matrix(pmc.output$ess, ncol = 1), 
            sprintf("results/drilling/pmc_%s_ess.txt", output.label),
            row.names = F, col.names = F)