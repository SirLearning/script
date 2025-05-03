library(mvtnorm); library(rgl)

par(mfrow = c(3, 2))
for (i in c(-0.75, 0, 0.75)) {
  rho = i; mu = c(0, 0); sig1 = sig2 = 1
  Sigma = matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), nrow = 2)
  X = seq(-3, 3, length = 41); Y = X
  f = function(X, Y) {
    XY = cbind(X, Y)
    dmvnorm(XY, mu, Sigma)
  }
  Z = outer(X, Y, f)
  persp(X, Y, Z, col = "lightgreen",
        theta = 30, phi = 20, r = 50, d = 0.1, expand = 0.5, 
        ltheta = 90, lphi = 180, shade = 0.75, ticktype = "detailed", nticks = 5,
        xlab = "X", ylab = "Y", zlab = "f(x,y)")
  # win.graph()
  image(X, Y, Z)
  contour(X, Y, Z, add = TRUE)
}

# 3D plot
persp3d(X, Y, Z, theta = 45, phi = 25, col = "red")
M = par3d("userMatrix"); bg3d("white")
play3d(par3dinterp(userMatrix = list(M, rotate3d(M, angle = pi,
                                                 x = 1, y = 0, z = 0))),
       duration = 10)