xdata <- data.frame(PH = related_data$PH,
                    HD = related_data$HD)
ydata <- data.frame(PH_Rep1 = plant_height$Rep1,
                    PH_Rep2 = plant_height$Rep2)
TSdiff.conf.ellipse = function(xdata, ydata, alpha) {
  if(ncol(xdata)!=2) stop("Only for bivariate normal")
  if(ncol(ydata)!=2) stop("Only for bivariate normal")
  n = nrow(xdata); xbar = colMeans(xdata); S1 = cov(xdata)
  m = nrow(ydata); ybar = colMeans(ydata); S2 = cov(ydata)
  S.pool = ((n-1)/(n+m-2))*S1+((m-1)/(n+m-2))*S2
  S = (1/n+1/m)*S.pool; es = eigen(S)
  e1 = es$vec %*% diag(sqrt(es$val))
  r1 = sqrt(qf(0.95, 2, n+m-3))*sqrt(2*(n+m-2)/(n+m-3))
  theta = seq(0, 2*pi, len = 250)
  v1 = cbind(r1*cos(theta), r1*sin(theta))
  pts = t(xbar-ybar-(e1%*%t(v1)))
  plot(pts, type = "l", col = 'red', lwd = 3,
       main = "Confidence Region for the difference of two mean vector", asp = 1)
  grid(lty = 1, equilogs = FALSE)
  th2 = c(0, pi/2, pi, 3*pi/2, 2*pi)
  v2 = cbind(r1*cos(th2), r1*sin(th2))
  pts2 = t(xbar-ybar-(e1%*%t(v2)))
  segments(pts2[3,1], pts2[3,2], pts2[1,1], pts2[1,2], lty = 4, lwd = 2)
  segments(pts2[2,1], pts2[2,2], pts2[4,1], pts2[4,2], lty = 4, lwd = 2)
}