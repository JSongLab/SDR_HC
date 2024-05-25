source("hcsdr.R")
#### Data generating
n = 400; p = 10; sd = 0.2;
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

# Run
## Modern - SDR
DCOV = sdr_via_dcov(X = x, Y = y)
HSIC = sdr_via_HSIC(X = x, Y = y)
BCOV = sdr_via_bcov(X = x, Y = y)

#### Measure of Performance : Modern - SDR
first_eigen(proj_matrix(beta) - proj_matrix(DCOV))
first_eigen(proj_matrix(beta) - proj_matrix(HSIC))
first_eigen(proj_matrix(beta) - proj_matrix(BCOV))


## HC - SDR
hcdr_DCOV = hcsdr_modern(x, y, "DCOV")
hcdr_HSIC = hcsdr_modern(x, y, "HSIC")
hcdr_BCOV = hcsdr_modern(x, y, "BCOV")

#### Measure of Performance : HC-SDR
diff_measure(beta,hcdr_DCOV)
diff_measure(beta,hcdr_HSIC)
diff_measure(beta,hcdr_BCOV)

#### Outcome direction
hcdr_DCOV
hcdr_HSIC
hcdr_BCOV
normalize(beta1)
