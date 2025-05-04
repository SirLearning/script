library(MASS); library(car)

set.seed(1); n = 10000
mu = c(0, 2); Sig = matrix(c(1, 0.7, 0.7, 1), 2, 2)
# method 1
biv1 = mvrnorm(n, mu, Sig); colnames(biv1) = c("X", "Y")
# method 2
sig.eigen = eigen(Sig)
Sigma = sig.eigen$vec%*%diag(sqrt(sig.eigen$val))%*%t(sig.eigen$vec)
biv2 = mu + t(matrix(rnorm(2*n), ncol = 2)%*%Sigma)
biv2 = t(biv2); colnames(biv2) = c("X", "Y")
# plot
par(mfrow = c(1, 2))
dataEllipse(biv1[, 1], biv1[, 2], levels = 0.95)
points(t(mu), col = "red", pch = 8)
dataEllipse(biv2[, 1], biv2[, 2], levels = 0.95)
points(t(mu), col = "red", pch = 8)
