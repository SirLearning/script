# Create data:
sampleNum = c(10, 50, 100, 200, 500, 1000, 2012)
# disc
disc = c(1.86, 1.44, 1.83, 2.06, 1.73, 1.48, 1.76)
blib = c(0.19, 0.14, 0.10, 0.07, 0.13, 0.25, 0.30)
scan = c(0.39, 0.53, 0.79, 1.12, 1.25, 0.98, 0.86)

par(pin = c(6, 6))

# # log mode
# plot(logRealTime~sampleNum, type="b", bty="l",
#      xlab="value of a", ylab="value of b", cex.lab = 1.5,
#      col=rgb(0.2,0.4,0.1,0.7), lwd=3, pch=17, ylim=c(3,10),
#      main = "Times", cex.main = 1.5
# )
# lines(logUserTime~sampleNum , col=rgb(0.8, 0.4, 0.1, 0.7) , lwd=10 , pch=19 , type="b" )
# lines(logSysTime~sampleNum , col=rgb(0.1, 0.6, 0.8, 0.7) , lwd=3 , pch=15 , type="b" )

# Make a basic graph
plot(disc~sampleNum, type="b", bty="l",
     xlab="sample number", ylab="instructions per CPU cycle (IPC)", cex.lab = 2,
     col=rgb(0.2,0.4,0.1,0.7), lwd=5, pch=17,
     ylim=c(min(c(disc, blib, scan)), max(c(disc, blib, scan))), cex.axis = 1.5,
     main = "IPC analysis", cex.main = 2
)
lines(blib~sampleNum , col=rgb(0.8, 0.4, 0.1, 0.7) , lwd=5 , pch=19 , type="b" )
lines(scan~sampleNum , col=rgb(0.1, 0.6, 0.8, 0.7) , lwd=5 , pch=15 , type="b" )

# Add a legend
legend("bottomleft",
       legend = c("disc", "blib", "scan"),
       col = c(rgb(0.2,0.4,0.1,0.7),
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1, 0.6, 0.8, 0.7)),
       pch = c(17,19,15),
       bty = "n",
       pt.cex = 2,
       cex = 1.5,
       text.col = "black",
       horiz = F ,
       inset = c(0.5, 0.8))