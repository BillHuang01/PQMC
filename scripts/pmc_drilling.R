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

# parameter setting
N <- 13
p <- 3
J <- 7
steps <- 7
sigma <- 0.2

# initialization
# library(mined)
# set.seed(950922)
# ini <- Lattice(N, p)
# ini[,1] <- 0.5 + 0.5 * ini[,1]
# ini[,2] <- 0.5 + 0.5 * ini[,2]
# ini[,3] <- 0.75 + 0.5 * ini[,3]
# pairs(ini)
# mean(dist(ini))
# write.table(ini, sprintf("data/drilling_ini_%d.txt", N), 
#             row.names = F, col.names = F, sep = "\t")


ini <- as.matrix(read.table(sprintf("data/drilling_ini_%d.txt", N)))

# Support Points
set.seed(950922)
pmc.sp <- pmc(ini, logf, J, steps, sigma, resample = "SP", output = 'raw',
              sigma.adapt = T, qmc = T, visualization = F)

write.table(pmc.sp$ini.all, "results/drilling/pmc_sp_ini.txt", row.names = F, col.names = F, sep = "\t")
write.table(pmc.sp$samp.all, "results/drilling/pmc_sp_samp.txt", row.names = F, col.names = F, sep = "\t")
write.table(matrix(pmc.sp$samp.all.logwts, ncol = 1), "results/drilling/pmc_sp_samplogwts.txt", row.names = F, col.names = F)
write.table(matrix(pmc.sp$ess, ncol = 1), "results/drilling/pmc_sp_ess.txt", row.names = F, col.names = F)

# Multinomial
set.seed(950922)
pmc.mn <- pmc(ini, logf, J, steps, sigma, resample = "Multinomial", output = 'raw',
              sigma.adapt = T, qmc = T, visualization = F)

write.table(pmc.mn$ini.all, "results/drilling/pmc_mn_ini.txt", row.names = F, col.names = F, sep = "\t")
write.table(pmc.mn$samp.all, "results/drilling/pmc_mn_samp.txt", row.names = F, col.names = F, sep = "\t")
write.table(matrix(pmc.mn$samp.all.logwts, ncol = 1), "results/drilling/pmc_mn_samplogwts.txt", row.names = F, col.names = F)
write.table(matrix(pmc.mn$ess, ncol = 1), "results/drilling/pmc_mn_ess.txt", row.names = F, col.names = F)

# Systematic
set.seed(950922)
pmc.ss <- pmc(ini, logf, J, steps, sigma, resample = "Systematic", output = 'raw',
              sigma.adapt = T, qmc = T, visualization = F)

write.table(pmc.ss$ini.all, "results/drilling/pmc_ss_ini.txt", row.names = F, col.names = F, sep = "\t")
write.table(pmc.ss$samp.all, "results/drilling/pmc_ss_samp.txt", row.names = F, col.names = F, sep = "\t")
write.table(matrix(pmc.ss$samp.all.logwts, ncol = 1), "results/drilling/pmc_ss_samplogwts.txt", row.names = F, col.names = F)
write.table(matrix(pmc.ss$ess, ncol = 1), "results/drilling/pmc_ss_ess.txt", row.names = F, col.names = F)