source("hcsdr.R")
## Data generating
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

## nonsparse
n = 400
scenario = 4
p <- 10
beta1 <- rep(1, 10) 
beta1 <- normalize(beta1)
beta2 <- c(1, 1, 1, -1, -1, -1, -1, 1, 1, -1)
beta2 <- normalize(beta2)
beta3 <- c(3, -1, 4, -2, -4, 5, 1, -3, -5, 2)
beta3 <- normalize(beta3)


x <- matrix(rnorm(n * p), n, p)

error <- rnorm(n, sd = sd)

beta <- switch(scenario, beta1, beta2, beta1, beta3)
y <- switch(scenario,
            (x %*% beta) + error,
            (x %*% beta)^2  + error,
            exp(x %*% beta) + error,
            5 * sin(x %*% beta) + error
)

# Run
# Run: Traditional SDR
dat <- data.frame(y = y,x = x)
SDRmethod = "DR"
dr_result = switch(SDRmethod,
                   "SIR" = dr::dr(y~., data = dat, method = "sir")$evec[,1],
                   "SAVE" = dr::dr(y~., data = dat, method = "save")$evec[,1],
                   "DR" = drr(x, y, ceiling(log(n)), 1),
                   "MAVE" = MAVE::mave(y~., data = dat, max.dim = 1)$dir[[1]],
                   "method error")

init <- normalize(dr_result)
# Run: HC SDR
hcdr = hcsdr(x,y, SDRmethod=SDRmethod)

# Measure of Performance
# Measure of Performance: HC-SDR
diff_measure(beta, hcdr)
# Measure of Performance: Traditional SDR
diff_measure(beta, init)

hcdr
beta1
