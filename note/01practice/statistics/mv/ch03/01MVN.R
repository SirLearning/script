library(mvtnorm)

# density function
dmvnorm(x, mean, sigma, log = FALSE)
# distribution function
pmvnorm(lower = Inf, upper = Inf, mean = rep(0, length(lower)),
        corr = NULL, sigma = NULL, algorithm = GenzBretz(), ...)
# quantile
qmvnorm(p, interval = NULL, tail = c("lower.tail", "upper.tail", "both.tails"),
        mean = 0, corr = NULL, sigma = NULL, algorithm = GenzBretz(), ...)
# random number of multivariable normal distribution
rmvnorm(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
        method = c("eigen", "svd", "chol"))