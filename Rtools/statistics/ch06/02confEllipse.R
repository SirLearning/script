xdata <- data.frame(PH = related_data$PH,
                    HD = related_data$HD)
alpha <- 0.95
pts <- function(xdata, alpha) {
  n = nrow(xdata); xbar = colMeans(xdata)
  S = cov(xdata); es = eigen(S)
  e1 = es$vec %*% diag(sqrt(es$val))
  r1 = sqrt(qf(alpha, 2, n-2)*sqrt(2*(n-1)/(n*(n-2))))
  theta = seq(0, 2*pi, len = 1000)
  v1 = cbind(r1*cos(theta), r1*sin(theta))
  pts = t(xbar - (e1%*%t(v1)))
  return(pts)
}

conf.elipse = function(xdata, alpha) {
  if (ncol(xdata) != 2) stop("Only for bivariate normal")
  pts_95 = pts(xdata, 0.95)
  pts_80 = pts(xdata, 0.8)
  plot(pts_95, type = "l", col = "red", lwd = 3,
       main = "Confidence Region for Bivariate Normal",
       xlab = "Plant Height", ylab = "Heading Date", asp = 1)
  lines(pts_80, type = "l", col = "yellow", lwd = 3)
  grid(lty = 1, equilogs = FALSE)
  segments(-0.2, xbar[2], xbar[1], xbar[2], lty = 2, lwd = 2)
  segments(xbar[1], 0, xbar[1], xbar[2], lty = 2, lwd = 2)
  th2 = c(0, pi/2, pi, 3*pi/2, 2*pi)
  v2 = cbind(r1*cos(th2), r1*sin(th2))
  pts2 = t(xbar - (e1%*%t(v2)))
  segments(pts2[3, 1], pts2[3, 2], pts2[1, 1], pts2[1, 2], lty = 4, lwd = 2)
  segments(pts2[2, 1], pts2[2, 2], pts2[4, 1], pts2[4, 2], lty = 4, lwd = 2)
}


