T2.test = function(X, mu0) {
  n = nrow(X); p = ncol(X)
  Xbar = colMeans(X); S = cov(X)
  T2 = n*t(Xbar-mu0)%*%solve(S)%*%(Xbar-mu0)
  T2.adj = (n-p)*T2/(p*(n-1))
  pval = 1 - pf(T2.adj, p, n-p)
  cat("Hotelling T-squared test for unknown covariance", fill = T)
  data.frame(T2 = T2, T2.adj = T2.adj, p.value = pval)
}

smp.PH <- data.frame(Rep1 = plant_height$Rep1,
                    Rep2 = plant_height$Rep2,
                    Rep3 = plant_height$Rep3,
                    Rep4 = plant_height$Rep4)
mu0 = c(100, 100, 100, 100)
T2.test(smp.PH, mu0)

# ICSNO method
library(ICSNP)
HotellingsT2(smp.PH, mu = mu0)

# Wilks
Y = factor(plant_height$country); X = as.matrix(plant_height[, 2:5])
fit = manova(X~Y)
summary(fit, test = "Wilks")
