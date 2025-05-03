library(ggplot2); library(ggiraphExtra)

ggRadar(data = plant_height, aes(group = country), alpha = 0) +
  theme(axis.text = element_text(size = 10), 
        legend.position = "right",
        legend.text = element_text(size = 10))

# cannot figure out why