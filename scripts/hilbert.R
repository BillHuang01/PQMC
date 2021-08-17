# hilbert curve package using bigz to handle large integers

# load bigz library
library(gmp)

# Decimal to Bianry
DecToBin <- function(decimal, width = 64){
  # input bigz object
  # output binary object
  if (!is.bigz(decimal)) stop("input must be a bigz object!")
  binary <- rep(0, width)
  n <- decimal
  i <- 1
  while (n > 0){
    r <- n %% 2
    binary[i] <- as.integer(r)
    n <- (n - r) / 2
    n <- as.bigz(n)
    i <- i + 1
  }
  binary <- rev(binary)
  return (binary)
}

# Binary to Decimal
BinToDec <- function(binary){
  # input binary object
  # output bigz object
  if (any(!(binary %in% c(0,1)))) stop("input must be a binary object")
  p <- length(binary)
  binary <- rev(binary)
  decimal <- as.bigz(0)
  for (i in 1:p){
    decimal <- decimal + as.bigz(2^(i-1)) * binary[i]
  }
  return (decimal)
}

hilbertcurve_coordinates_from_distance <- function(d, dim, order){
  #################################################
  # variable 
  # d: distance along hilbert curve (double, 0 - 1)
  # dim: dimensions (integer)
  # order: iterations of the hilbert curve (integer)
  #################################################
  
  # condition check
  if (d < 0 |d > 1) stop("distance must be normalied to 0 and 1!")
  if (dim < 1) stop("dimension must be positive!")
  if (order < 1) stop("order must be positive!")
  
  # convert to coordinate
  d <- d * as.bigq(as.bigz(2^(order*dim)))
  d <- as.bigz(d)
  d.bit <- DecToBin(d, order*dim)
  x <- rep(0,dim)
  for (i in 1:dim){
    x[i] <- as.integer(BinToDec(d.bit[seq(from=i,by=dim,length.out=order)]))
  }
  
  # Gray decode
  t <- bitwShiftR(x[dim],1) # >> operator in python
  for (i in dim:2) x[i] <- bitwXor(x[i], x[(i-1)]) # ^ operator in python
  x[1] <- bitwXor(x[1], t)
  
  # Undo excess work
  Z = bitwShiftL(2,(order-1))
  Q = 2
  while (Q != Z){
    P <- Q - 1
    for (i in dim:1){
      if (bitwAnd(x[i],Q)){
        x[1] <- bitwXor(x[1], P) # invert
      } else{
        t <- bitwAnd(bitwXor(x[1], x[i]), P)
        x[1] <- bitwXor(x[1], t) # exchange
        x[i] <- bitwXor(x[i], t) # exchange
      }
    }
    Q = bitwShiftL(Q,1)
  }
  
  return (x)
  
}

hilbertcurve_distance_from_coordinates <- function(x, dim, order){
  #################################################
  # variable 
  # x: coordinate (integer)
  # dim: dimensions (integer)
  # order: iterations of the hilbert curve (integer)
  #################################################
  
  # condition check
  if (dim < 1) stop("dimension must be positive!")
  if (order < 1) stop("order must be positive!")
  if (length(x)!=dim) stop(sprintf("x does not have dimension %d!", dim))
  if (any(x<0)) stop("invalid coordinate input, one of more dimensions 
                      have value less than 0!")
  if (any(x>(2^order-1))) stop("invalid coordinate input, one or mode dimensions 
                                have value greater than 2^order - 1!")
  
  # inverse undo excess work
  M <- bitwShiftL(1, order-1)
  Q <- M
  while (Q > 1){
    P <- Q - 1
    for (i in 1:dim){
      if (bitwAnd(x[i],Q)){
        x[1] <- bitwXor(x[1], P)
      } else {
        t <- bitwAnd(bitwXor(x[1],x[i]), P)
        x[1] <- bitwXor(x[1], t)
        x[i] <- bitwXor(x[i], t)
      }
    }
    Q <- bitwShiftR(Q, 1)
  }
  
  # Gray encode
  for (i in 2:dim) x[i] = bitwXor(x[i], x[(i-1)])
  t <- 0
  Q <- M
  while (Q > 1){
    if (bitwAnd(x[dim],Q)) t <- bitwXor(t,Q-1)
    Q <- bitwShiftR(Q, 1)
  }
  for (i in 1:dim) x[i] <- bitwXor(x[i], t)
  
  # convert it to distance
  x.bit <- rep(0, order*dim)
  for (i in 1:dim){ 
    x.bit[seq(from=i,by=dim,length.out=order)] <- DecToBin(as.bigz(x[i]), order)
  }
  d <- BinToDec(x.bit)
  d <- as.numeric(d / as.bigz(2^(order*dim)))

  return (d)
  
}
