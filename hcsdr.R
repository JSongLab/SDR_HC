# HCR subcodes
require(HellCor)
require(SphereOptimize)
require(dr)
require(MAVE)
require(beepr)
require(doParallel)
require(energy)

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

iht <- function(x, y, r){
  z <- standmat(x)
  szy <- cov(z, y)
  szz <- var(z)
  p <- dim(z)[2]
  imat <- szy
  for(i in 1:(p - 1)){
    imat <- cbind(imat, matpower(szz,i) %*% szy)
  }
  return(eigen(imat %*% t(imat))$vectors[, 1:r])
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