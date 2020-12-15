Population Quasi-Monte Carlo
================
Chaofan (Bill) Huang

Let us first load the required libraries and scripts.

``` r
set.seed(090820)
source("scripts/lib.R")
library(mvtnorm)
library(randtoolbox)
```

Now define the log density for the two dimensional mixture of five
normals.

``` r
logmixture <- function(x){
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
```

Let us draw the density contour of the mixture distribution.

``` r
x1 <- x2 <- seq(0, 1, length.out = 101)
x.grid <- expand.grid(x1, x2)
density <- matrix(exp(apply(x.grid, 1, logmixture)), 101, 101)
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15)
```

![](README_files/figure-gfm/contour-1.png)<!-- -->

## Resampling Method Comparsion

Consider resmaple n = 25 points from 2500 Sobol points over the unit
square as importance samples for the mixture of normals using
Multinomial, Residual, Stratified, Systematic, and Support Points.

``` r
n <- 25
N <- 2500
layout(matrix(c(1:6), nrow = 2, byrow = T))
# Sobol points
samp <- sobol(N, 2)
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15, main = "Sobol")
points(samp, pch = 18, cex = 0.5, col = "green")
# compute the weight
samp.logwts <- apply(samp, 1, logmixture)
samp.wts <- exp(samp.logwts - max(samp.logwts))
samp.wts <- samp.wts / sum(samp.wts)
# Multinomial resampling
mn.samp <- samp[sample(1:N, n, replace=T, prob=samp.wts),]
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15, main = "Multinomial")
points(mn.samp, pch = 16, cex = 1, col = "red")
# Residual resampling
rs.samp <- samp[rs.sample(1:N, n, prob=samp.wts),]
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15, main = "Residual")
points(mn.samp, pch = 16, cex = 1, col = "red")
# Stratified resampling
st.samp <- samp[st.sample(1:N, n, prob=samp.wts),]
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15, main = "Stratified")
points(st.samp, pch = 16, cex = 1, col = "red")
# Systematic resampling
ss.samp <- samp[ss.sample(1:N, n, prob=samp.wts),]
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15, main = "Systematic")
points(ss.samp, pch = 16, cex = 1, col = "red")
# Importance Support Points resampling
sp.samp <- samp[sp.sample(samp, n, prob=samp.wts),]
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15, main = "ISP")
points(sp.samp, pch = 16, cex = 1, col = "red")
```

![](README_files/figure-gfm/resample-1.png)<!-- -->

## PQMC vs. PMC on Mixture of Normals

Let us now run PMC and PQMC on the mixture of five normals with N = 50,
J = 20, and T = 6. The initial centers are the 50 Sobol points over the
unit square. Adaptation for covariance is applied. Here follows the
parameter setting.

``` r
expectation <- 0.5 + c(1.6,1.4) / 40 # E[X]
Z <- 1 # normalizing constant
p <- 2 # dimensions
N <- 50 # number of propopsals
J <- 20 # number of samples drawn from each proposal
steps <- 6 # number of PMC iterations
# initilization
ini <- sobol(N, p) # initial centers
ini.logq <- NULL
# variance
sigma <- 0.2 # initial covariance
sigma.adapt <- T # adaptation for covariance
```

### PMC (Multinomial)

``` r
layout(matrix(c(1:6), nrow = 2, byrow = T))
pmc.mn <- pmc(logmixture, N, J, steps, ini, ini.logq,
              sampling = "random", resampling = "multinomial",
              sigma = sigma, sigma.adapt = sigma.adapt, visualization = T)
```

![](README_files/figure-gfm/pmc-mixture-mn-1.png)<!-- -->

``` r
# log MSE of E[X] by standard PMC estimator
log(mean((pmc.mn$m.std - expectation)^2))
```

    ## [1] -8.466343

``` r
# log MSE of E[X] by weighted PMC estimator
log(mean((pmc.mn$m.wts - expectation)^2))
```

    ## [1] -9.724102

``` r
# log MSE of Z by standard PMC estimator
log((pmc.mn$z.std - Z)^2)
```

    ## [1] -6.652152

``` r
# log MSE of Z by weighted PMC estimator
log((pmc.mn$z.wts - Z)^2)
```

    ## [1] -7.833203

### PMC (Systematic)

``` r
layout(matrix(c(1:6), nrow = 2, byrow = T))
pmc.ss <- pmc(logmixture, N, J, steps, ini, ini.logq,
              sampling = "random", resampling = "systematic",
              sigma = sigma, sigma.adapt = sigma.adapt, visualization = T)
```

![](README_files/figure-gfm/pmc-mixture-ss-1.png)<!-- -->

``` r
# log MSE of E[X] by standard PMC estimator
log(mean((pmc.ss$m.std - expectation)^2))
```

    ## [1] -10.39755

``` r
# log MSE of E[X] by weighted PMC estimator
log(mean((pmc.ss$m.wts - expectation)^2))
```

    ## [1] -13.00525

``` r
# log MSE of Z by standard PMC estimator
log((pmc.ss$z.std - Z)^2)
```

    ## [1] -7.619664

``` r
# log MSE of Z by weighted PMC estimator
log((pmc.ss$z.wts - Z)^2)
```

    ## [1] -7.761437

### PQMC (Importance Support Points)

``` r
layout(matrix(c(1:6), nrow = 2, byrow = T))
pqmc.sp <- pmc(logmixture, N, J, steps, ini, ini.logq,
               sampling = "qmc", resampling = "sp",
               sigma = sigma, sigma.adapt = sigma.adapt, visualization = T)
```

![](README_files/figure-gfm/pqmc-mixture-sp-1.png)<!-- -->

``` r
# log MSE of E[X] by standard PMC estimator
log(mean((pqmc.sp$m.std - expectation)^2))
```

    ## [1] -9.772567

``` r
# log MSE of E[X] by weighted PMC estimator
log(mean((pqmc.sp$m.wts - expectation)^2))
```

    ## [1] -13.40994

``` r
# log MSE of Z by standard PMC estimator
log((pqmc.sp$z.std - Z)^2)
```

    ## [1] -7.680051

``` r
# log MSE of Z by weighted PMC estimator
log((pqmc.sp$z.wts - Z)^2)
```

    ## [1] -10.24683

## PQMC vs. PMC on Banana Shape Distribution

Here is the contour of the banana shape distribution.

``` r
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
x1 <- x2 <- seq(0, 1, length.out = 101)
x.grid <- expand.grid(x1, x2)
density <- matrix(exp(apply(x.grid, 1, logbanana)), 101, 101)
contour.default(x = x1, y = x2, z = density, drawlabels = F, nlevels = 15)
```

![](README_files/figure-gfm/banana-1.png)<!-- -->

Let us now run PMC and PQMC on the mixture of five normals with N = 50,
J = 20, and T = 6. The initial centers are the 50 Sobol points over the
unit square. Adaptation for covariance is applied. Here follows the
parameter setting.

``` r
expectation <- c(0.5,0.7162834) # E[X]
Z <- 0.0223886 # normalizing constant
# Above computed from grid approximation
p <- 2 # dimensions
N <- 50 # number of propopsals
J <- 20 # number of samples drawn from each proposal
steps <- 6 # number of PMC iterations
# initilization
ini <- sobol(N, p) # initial centers
ini.logq <- NULL
# variance
sigma <- 0.2 # initial covariance
sigma.adapt <- T # adaptation for covariance
```

### PMC (Multinomial)

``` r
layout(matrix(c(1:6), nrow = 2, byrow = T))
pmc.mn <- pmc(logbanana, N, J, steps, ini, ini.logq,
              sampling = "random", resampling = "multinomial",
              sigma = sigma, sigma.adapt = sigma.adapt, visualization = T)
```

![](README_files/figure-gfm/pmc-banana-mn-1.png)<!-- -->

``` r
# log MSE of E[X] by standard PMC estimator
log(mean((pmc.mn$m.std - expectation)^2))
```

    ## [1] -6.388404

``` r
# log MSE of E[X] by weighted PMC estimator
log(mean((pmc.mn$m.wts - expectation)^2))
```

    ## [1] -8.051862

``` r
# log MSE of Z by standard PMC estimator
log((pmc.mn$z.std - Z)^2)
```

    ## [1] -15.99496

``` r
# log MSE of Z by weighted PMC estimator
log((pmc.mn$z.wts - Z)^2)
```

    ## [1] -13.40895

### PMC (Systematic)

``` r
layout(matrix(c(1:6), nrow = 2, byrow = T))
pmc.ss <- pmc(logbanana, N, J, steps, ini, ini.logq,
              sampling = "random", resampling = "systematic",
              sigma = sigma, sigma.adapt = sigma.adapt, visualization = T)
```

![](README_files/figure-gfm/pmc-banana-ss-1.png)<!-- -->

``` r
# log MSE of E[X] by standard PMC estimator
log(mean((pmc.ss$m.std - expectation)^2))
```

    ## [1] -10.54613

``` r
# log MSE of E[X] by weighted PMC estimator
log(mean((pmc.ss$m.wts - expectation)^2))
```

    ## [1] -11.10372

``` r
# log MSE of Z by standard PMC estimator
log((pmc.ss$z.std - Z)^2)
```

    ## [1] -17.06236

``` r
# log MSE of Z by weighted PMC estimator
log((pmc.ss$z.wts - Z)^2)
```

    ## [1] -20.05266

### PQMC (Importance Support Points)

``` r
layout(matrix(c(1:6), nrow = 2, byrow = T))
pqmc.sp <- pmc(logbanana, N, J, steps, ini, ini.logq,
               sampling = "qmc", resampling = "sp",
               sigma = sigma, sigma.adapt = sigma.adapt, visualization = T)
```

![](README_files/figure-gfm/pqmc-banana-sp-1.png)<!-- -->

``` r
# log MSE of E[X] by standard PMC estimator
log(mean((pqmc.sp$m.std - expectation)^2))
```

    ## [1] -11.42542

``` r
# log MSE of E[X] by weighted PMC estimator
log(mean((pqmc.sp$m.wts - expectation)^2))
```

    ## [1] -12.03573

``` r
# log MSE of Z by standard PMC estimator
log((pqmc.sp$z.std - Z)^2)
```

    ## [1] -15.31548

``` r
# log MSE of Z by weighted PMC estimator
log((pqmc.sp$z.wts - Z)^2)
```

    ## [1] -15.69704
