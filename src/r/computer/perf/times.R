# Create data:
sampleNum = c(10, 50, 100, 200, 500, 1000, 2012)
# disc
realTime = c(22.412912940, 61.987440044, 67.148435684, 98.954823747, 275.691697084, 647.531322019, 1080.447711542)
sysTime = c(33.384027000, 255.297659000, 342.131588000, 411.730050000, 971.025948000, 3537.768358000, 2955.494378000)
userTime = c(176.303062000, 1381.087839000, 2152.393975000, 4085.091478000, 12836.716418000, 29399.571239000, 56254.484098000)

# blib
# realTime = c(3.749421021, 5.847891564, 8.030426103, 12.581682366, 17.041002394, 42.514106507, 79.451093958)
# sysTime = c(2.192140000, 11.580206000, 23.451592000, 23.934930000, 28.019004000, 164.035906000, 286.456827000)
# userTime = c(51.798211000, 104.373806000, 174.630179000, 352.061580000, 275.715227000, 1799.116253000, 2469.696131000)

# scan
# realTime = c(6.948137147, 18.235428376, 29.329003329, 59.369517816, 183.530870764, 580.178783118, 1852.264429938)
# sysTime = c(3.334818000, 41.394046000, 53.838342000, 196.938010000, 959.080038000, 3910.914478000, 10667.065494000)
# userTime = c(129.853238000, 762.789135000, 1257.880222000, 2079.958141000, 5811.457374000, 20140.568745000, 79400.466260000)

# log mode
# logRealTime = log(realTime)
# logUserTime = log(userTime)
# logSysTime = log(sysTime)

par(pin = c(4, 6))
linters: with_defaults(line_length_linter = NULL)
# # log mode
# plot(logRealTime~sampleNum, type="b", bty="l",
#      xlab="value of a", ylab="value of b", cex.lab = 1.5,
#      col=rgb(0.2,0.4,0.1,0.7), lwd=3, pch=17, ylim=c(3,10),
#      main = "Times", cex.main = 1.5
# )
# lines(logUserTime~sampleNum , col=rgb(0.8, 0.4, 0.1, 0.7) , lwd=10 , pch=19 , type="b" )
# lines(logSysTime~sampleNum , col=rgb(0.1, 0.6, 0.8, 0.7) , lwd=3 , pch=15 , type="b" )

# Make a basic graph
plot(realTime~sampleNum, type="b", bty="l",
     xlab="sample number", ylab="time cost (seconds)", cex.lab = 2,
     col=rgb(0.2,0.4,0.1,0.7), lwd=5, pch=17,
     ylim=c(min(c(realTime, userTime, sysTime)), max(c(realTime, userTime, sysTime))), cex.axis = 1.5,
     main = "disc run times", cex.main = 2
)
lines(userTime~sampleNum , col=rgb(0.8, 0.4, 0.1, 0.7) , lwd=5 , pch=19 , type="b" )
lines(sysTime~sampleNum , col=rgb(0.1, 0.6, 0.8, 0.7) , lwd=5 , pch=15 , type="b" )

# Add a legend
legend("bottomleft",
       legend = c("real time", "user time", "system time"),
       col = c(rgb(0.2,0.4,0.1,0.7),
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1, 0.6, 0.8, 0.7)),
       pch = c(17,19,15),
       bty = "n",
       pt.cex = 2,
       cex = 1.5,
       text.col = "black",
       horiz = F ,
       inset = c(0.1, 0.8))