library(MASS); library(ade4); library(scatterplot3d)

# color
rb = rainbow(length(unique(plant_height$country)))
rb = adjustcolor(rb, alpha.f = 0.4)


# panel plot
panel1 = function(X, Y) {
  XY = cbind.data.frame(X, Y)
  s.class(XY, as.factor(plant_height$country), 
          include.ori = F, add.p = T, clab = 0, col = rb, cpoi = 2, csta = 0)
}
pairs(plant_height[, 2:5], panel = panel1)

# 3D scatter
par(mfrow = c(2, 2)); mar0 = c(3, 3, 1, 3)
scatterplot3d(plant_height[, 2], plant_height[, 3], plant_height[, 4],
              mar = mar0, color = rb[as.factor(plant_height$country)], pch = 19,
              xlab = "Rep1", ylab = "Rep2", zlab = "Rep3")
scatterplot3d(plant_height[, 3], plant_height[, 4], plant_height[, 5],
              mar = mar0, color = rb[as.factor(plant_height$country)], pch = 19,
              xlab = "Rep2", ylab = "Rep3", zlab = "Rep4")
scatterplot3d(plant_height[, 2], plant_height[, 4], plant_height[, 5],
              mar = mar0, color = rb[as.factor(plant_height$country)], pch = 19,
              xlab = "Rep1", ylab = "Rep3", zlab = "Rep4")
scatterplot3d(plant_height[, 2], plant_height[, 3], plant_height[, 5],
              mar = mar0, color = rb[as.factor(plant_height$country)], pch = 19,
              xlab = "Rep1", ylab = "Rep2", zlab = "Rep4")


