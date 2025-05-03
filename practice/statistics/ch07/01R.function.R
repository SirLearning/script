lm.reg = lm(PH ~ HD, data = related_data)
summary(lm.reg)
confint(lm.reg)
library(GGally)
ggcoef(lm.reg, exclude_intercept = T, vline_color = "red",
       errorbar_color = "blue", errorbar_height = 0.1)+
  theme_bw()
y.res = residuals(lm.reg)
shapiro.test(y.res)
y.fit = predict(lm.reg)
plot(y.res~y.fit)
y.rst = rstandard(lm.reg)
plot(y.rst~y.fit)

plot_residuals = function(x) {
  plot(x, which = 1:6,
       caption = c("Residuals vs Fitted", "Normal Q-Q plot",
                   "Scale-Location plot", "Cook's distance plot"),
       panel = points, sub.caption = deparse(x$call), main = "",
       ask = prod(par("mfcol"))<length(which)&&dev.interactive(),
       id.n = 3, label.id = names(residuals(x)), cex.id = 0.75)
}
