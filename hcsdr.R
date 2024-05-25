# HCR subcodes
require(HellCor)
require(SphereOptimize)
require(dr)
require(MAVE)
require(doParallel)
require(Ball)
require(dHSIC)
require(energy)
require(Rdonlp2)


## Hellinger Correlation 
HCXY <- function(beta, X, Y){
  
  ### ERROR
  if(sum(is.nan(beta)) > 0) return(10000)
  
  xb  <- X %*% beta 
  z   <- cbind(xb, Y)
  n   <- dim(z)[1]
  rz  <- apply(z, 2, rank) / (n + 1) 
  u   <- get.knn(rz, 1)$nn.dist
  res <- mean(u) * sqrt(n - 1) * 2  
  return(res)
}

## projection matrix
proj_matrix <- function(mat){
  mat %*% solve(t(mat) %*% mat) %*% t(mat)
}

## largest eigenvalue
first_eigen <- function(mat){
  mat <- (mat + t(mat)) / 2
  eigval <- eigen(mat, symmetric = TRUE)$values
  max(eigval)
}


normalize <- function(vec){
  vec / sqrt(sum(vec^2))
}

matpower <- function(mat, alpha){
  mat <- round((mat + t(mat))/2, 7)
  tmp <- eigen(mat)
  tmp$vectors %*% diag((tmp$values)^alpha) %*% t(tmp$vectors)
}

discretize <- function(y, h){
  n <- length(y)
  m <- floor(n / h)
  y <- y + 0.00001 * mean(y) * rnorm(n)
  yord <- y[order(y)]
  divpt <- numeric()
  for(i in 1:(h - 1)){
    divpt = c( divpt, yord[i * m + 1])
  }
  y1 <- rep(0, n)
  y1[y < divpt[1]] <- 1
  y1[y >= divpt[h - 1]] <- h
  for(i in 2:(h - 1)){
    y1[(y >= divpt[i - 1]) & (y < divpt[i])] <- i
  }
  return(y1)
}

standmat <- function(mat){
  mu <- apply(mat, 2, mean)
  sig <- var(mat)
  signrt <- matpower(sig, -1/2)
  t(t(mat) - mu) %*% signrt
}

drr <- function(x, y, h, r){
  p <- ncol(x)
  n <- nrow(x)
  signrt <- matpower(var(x), -1/2)
  xc <- t(t(x) - apply(x, 2, mean))
  xst <- xc %*% signrt
  ydis <- discretize(y, h)
  yless <- ydis
  ylabel <- numeric()
  for(i in 1:n){
    if(var(yless) != 0){ 
      ylabel = c(ylabel, yless[1])
      yless <- yless[yless!=yless[1]]
    }
  }
  ylabel <- c(ylabel, yless[1])
  prob = numeric()
  for(i in 1:h){
    prob <- c(prob, length(ydis[ydis == ylabel[i]])/n)
  }
  vxy <- array(0, c(p, p, h))
  exy <- numeric()
  for(i in 1:h){
    vxy[, , i] <- var(xst[ydis == ylabel[i], ])
    exy <- rbind(exy, apply(xst[ydis == ylabel[i], ], 2, mean))
  }
  mat1 <- matrix(0, p, p)
  mat2 <- matrix(0, p, p)
  for(i in 1:h){
    mat1 <- mat1 + prob[i] * (vxy[, , i] + exy[i, ] %*% t(exy[i,])) %*% (vxy[, , i] + exy[i, ] %*% t(exy[i, ]))
    mat2 <- mat2 + prob[i] * exy[i, ] %*% t(exy[i, ])
  }
  out <- 2 * mat1 + 2 * mat2 %*% mat2 + 2 * sum(diag(mat2)) * mat2 - 2 * diag(p)
  return(signrt %*% eigen(out)$vectors[, 1:r])
}

diff_measure <- function(vec1, vec2){
  return(first_eigen(proj_matrix(vec1) - proj_matrix(vec2)))
}

## optimize
sphoptim <- function (par, fn, ...){
  ## Check initial value is on a unit sphere
  if (round(sum(par^2), 10) != 1) {
    stop("Initial value is not on a unit sphere.")
  }
  
  ## define obj function on sphere
  temp_fn <- function(t, ...) {
    s <- from.Sphere(t)
    res <- fn(s, ...)
    return(res)
  }
  
  ## change coordinate and optimize
  theta <- to.Sphere(par)
  k <- stats::optim(theta, temp_fn,...)
  
  ## back to original coordinate
  par <- from.Sphere(k$par)
  value <- k$value  
  count <- k$counts
  return(list(par = par, value = value, count = count))
}

hcsdr = function(x, y, SDRmethod="DR", maxit = 10000){
  f = function(beta) HCXY(beta,X=x,Y=y)
  drout <- switch(SDRmethod,
                  "SIR" = dr::dr(y~., data = dat, method = "sir")$evec[,1],
                  "SAVE" = dr::dr(y~., data = dat, method = "save")$evec[,1],
                  "DR" = drr(x, y, ceiling(log(n)), 1),
                  "MAVE" = MAVE::mave(y~., data = dat, max.dim = 1)$dir[[1]],
                  "method error")
  init <- normalize(drout)
  
  hcdr2 = SphereOptimize(init, f, method="SANN", control = list(maxit = maxit*2))$par
  hcdr = sphoptim(hcdr2, HCXY, X = x, Y = y,control = list(maxit = 10000))$par
  return(hcdr)
}

##########################
####  Modern SDR part ####
##########################

nlcon = function(vec) as.numeric(vec%*%vec)
dnlcon = function(vec) 2 * vec
attr(nlcon, "gr") =  dnlcon
nlin.l = nlin.u = 1

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
                  "DCOV" = sdr_via_dcov(X = x, Y = y),
                  "HSIC" = sdr_via_HSIC(X = x, Y = y),
                  "BCOV" = sdr_via_bcov(X = x, Y = y),
                  "method error")
  init <- normalize(drout)
  
  #### HCOR part
  hcdr2 = SphereOptimize(init, f, method="SANN", control = list(maxit = maxit*2))$par
  hcdr = sphoptim(hcdr2, HCXY, X = x, Y = y,control = list(maxit = 10000))$par
  return(hcdr)
}

