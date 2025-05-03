library(ggplot2); library(GGally)

ggparcoord(Watkins_phenotype, columns = 10:13, groupColumn = 3,
           scale = "globalminmax", showPoints = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  labs(x = "phenotype", y = "value")

ggparcoord(plant_height, columns = 2:5, groupColumn = 1,
           scale = "globalminmax", showPoints = TRUE) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  labs(x = "phenotype", y = "value")
