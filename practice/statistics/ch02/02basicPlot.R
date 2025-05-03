# par(mai = c(0.7,0.7,0.4,0.4), cex = 0.8)
data$PH_M_cm_Rep1 <- as.numeric(as.character(data$PH_M_cm_Rep1))
data$PH_M_cm_Rep2 <- as.numeric(as.character(data$PH_M_cm_Rep2))
plot(data$PH_M_cm_Rep1, data$PH_M_cm_Rep2, xlab = "hight 1", ylab = "hight 2",
     pch = 19, cex = 1.3, col = 'red')
grid(col = "grey60"); axis(side = 4, lty = 1)
points(mean(data$PH_M_cm_Rep1), mean(data$PH_M_cm_Rep2), pch = 19, cex = 4, col = 'black')
abline(v = mean(data$PH_M_cm_Rep1), h = mean(data$PH_M_cm_Rep2), lty = 2, col = "gray30")
abline(lm(data$PH_M_cm_Rep1~data$PH_M_cm_Rep2), lwd = 2, col = 'blue')
fit = lm(data$PH_M_cm_Rep1~data$PH_M_cm_Rep2)
data$predeited = predict(fit); data$residuals = residuals(fit)
segments(data$PH_M_cm_Rep1, data$PH_M_cm_Rep2, data$PH_M_cm_Rep1, data$predeited)
legend("topleft", legend = "regression line", lty = c(1,6),
       col = 'blue', cex = 0.95, fill = 'blue', box.col = "grey60",
       ncol = 1, inset = 0.01, x.intersp = 0.3)
box(col = 1, lwd = 2)
