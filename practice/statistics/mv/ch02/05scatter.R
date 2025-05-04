library(ggpubr); library(GGally); library(dplyr)

ggscatter(data = plant_height, x = "Rep1", y = "Rep4",
          add = "reg.line", conf.int = TRUE) +
  stat_regline_equation(label.x = 51, label.y = 130, size = 4) +
  stat_cor(label.x = 51, label.y = 135, size = 4) +
  theme_bw()

ggpairs(data = plant_height, columns = 2:5)
