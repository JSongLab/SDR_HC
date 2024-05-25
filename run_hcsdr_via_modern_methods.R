require(Ball)
require(dHSIC)
require(Rdonlp2)

source("hcsdr.R")
#### Data generating
n = 400; p = 10; sd=0.2;
x <-matrix(rnorm(n * p), n, p)
error <- rnorm(n, sd = sd)

scenario = 1
beta1 <- c(rep(1, 5), rep(0, 5))
beta2 <- c(1, -1, rep(0, 8))
beta3 <- c(1, 1, rep(0, 8))
beta <- switch(scenario, beta1, beta2, beta1, beta3)

y <- switch(scenario,
            (x %*% beta) + error,
            (x %*% beta)^2  + error,
            exp(x %*% beta) + error,
            5 * sin(x %*% beta) + error
)

#### Initial value for optimization for each dependency measure
init_value = function(x, y, SDRmethod = "DCOV", n = 200){
  p = dim(x)[2]
  X <<- x; Y <<- y
  value = 0
  drout = rep(0, p)
  for(i in 1:n){
    temp_dr = normalize(rnorm(p))
    temp_value = switch(SDRmethod,
                   "DCOV" = dcov_ftn(temp_dr),
                   "HSIC" = HSIC_ftn(temp_dr),
                   "BCOV" = bcov_ftn(temp_dr),
                   "method error")
    if(temp_value < value){
      value = temp_value
      drout = temp_dr
    }
  }
  return(drout)
}

#### DCOV part
dcov_ftn = function(init){
  xb = as.vector(X%*%init)
  y = as.vector(Y)
  return(-dcov(xb, y))
}

sdr_via_dcov = function(init = NULL, X, Y){
  if (is.null(init)){
    init = init_value(X, Y, "DCOV")
  }
  X <<- X; Y <<- Y
  p = dim(X)[2]
  par.l = rep(-1, p)
  par.u = rep(1, p)
  res = donlp2(init, dcov_ftn, par.l = par.l, par.u = par.u, nlin =list(nlcon), nlin.u = nlin.u, nlin.l = nlin.l)
  return(res$par)
}

#### HSIC part
HSIC_ftn = function(init){
  xb = as.vector(X%*%init)
  y = as.vector(Y)
  return(-dhsic(xb, y)$dHSIC)
}

sdr_via_HSIC = function(init = NULL, X, Y){
  if (is.null(init)){
    init = init_value(X, Y, "HSIC")
  }
  X <<- X; Y <<- Y
  p = dim(X)[2]
  par.l = rep(-1, p)
  par.u = rep(1, p)
  res = donlp2(init, HSIC_ftn, par.l = par.l, par.u = par.u, nlin =list(nlcon), nlin.u = nlin.u, nlin.l = nlin.l)
  return(res$par)
}

#### BCOV part
bcov_ftn = function(init){
  xb = as.vector(X%*%init)
  y  = as.vector(Y)
  return(-bcov(xb, y)) 
}


sdr_via_bcov = function(init = NULL, X, Y){
  if (is.null(init)){
    init = init_value(X, Y, "BCOV")
  }
  X <<- X; Y <<- Y
  p = dim(X)[2]
  par.l = rep(-1, p)
  par.u = rep(1, p)
  res = donlp2(init, bcov_ftn, par.l = par.l, par.u = par.u, nlin=list(nlcon), nlin.u=nlin.u, nlin.l=nlin.l)
  return(res$par)
}

hcsdr_modern = function(x, y, SDRmethod="DCOV", maxit = 10000){
  
  f = function(beta) HCXY(beta,X = x,Y = y)
  p = dim(x)[2]
  
  #### Optimization Constraint for Rdonlp2 package
  nlcon = function(vec) as.numeric(vec%*%vec)
  dnlcon = function(vec) 2 * vec
  attr(nlcon, "gr") =  dnlcon
  nlin.l = nlin.u = 1
  par.l = rep(-1, p)
  par.u = rep(1, p)
  
  drout <- switch(SDRmethod,
                  "DCOV" = sdr_via_dcov(init, X = x, Y = y),
                  "HSIC" = sdr_via_HSIC(init, X = x, Y = y),
                  "BCOV" = sdr_via_bcov(init, X = x, Y = y),
                  "method error")
  init <- normalize(drout)
  
  #### HCOR part
  hcdr2 = SphereOptimize(init, f, method="SANN", control = list(maxit = maxit*2))$par
  hcdr = sphoptim(hcdr2, HCXY, X = x, Y = y,control = list(maxit = 10000))$par
  return(hcdr)
}

hcdr_DCOV = hcsdr_modern(x, y, "DCOV")
hcdr_HSIC = hcsdr_modern(x, y, "HSIC")
hcdr_BCOV = hcsdr_modern(x, y, "BCOV")

#### Measure of Performance : HC-SDR
first_eigen(proj_matrix(beta1) - proj_matrix(hcdr_DCOV))
first_eigen(proj_matrix(beta1) - proj_matrix(hcdr_HSIC))
first_eigen(proj_matrix(beta1) - proj_matrix(hcdr_BCOV))

#### Outcome direction
hcdr_DCOV
hcdr_HSIC
hcdr_BCOV
normalize(beta1)
