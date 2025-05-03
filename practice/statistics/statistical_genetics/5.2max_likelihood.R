x=3
n=15
theta=seq(.1, .9, by=.01)
logl=c()
num=length(theta)
for (i in 1:num) {
  logl[i]=-2*log(dbinom(x, n, theta[i]))
}
plot(theta, logl, type="l",
     xlab="Theta", ylab="-2loglikelihood",
     font.lab=2, cex.lab=1.3, las=1, lwd=2,
     cex.axis=1.2)
minl=min(logl)
est=theta[logl==minl]
abline(h=minl, v=est)
grid(col = "darkgray")
print(c(minl, est))
