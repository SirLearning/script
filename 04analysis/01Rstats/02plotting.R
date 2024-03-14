library(ggsci)
# length of the genome
ggplot(data = allTE_summ, aes(x = Classification, y = size)) +
  geom_col(aes(fill = Clasification, group = Species),
           width = 0.7, color = "black", position = "dodge") +
  scale_fill_igv() +
  labs(x = NULL, color = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
# proportion of the genome
ggplot(data = allTE_summ, aes(x = Classification, y = percent)) +
  geom_col(aes(fill = Classification, group = Species),
           width = 0.7, color = "black", position = "dodge") +
  scale_fill_igv() +
  labs(x = NULL, color = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
# number of TEs
ggplot(data = intactTE_summ, aes(x = Classification, y = count)) +
  geom_col(aes(fill = Classification, group = Species),
           width = 0.7, color = "black", position = "dodge") +
  scale_fill_igv() +
  labs(x = NULL, color = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")