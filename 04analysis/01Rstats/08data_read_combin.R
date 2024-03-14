library(tidyverse)
Aly.chr_sw_dencity = bind_rows(
  read_table(
    file = "34.LTR_chr//Aly/chr_sw.gene.density.txt") %>%
    mutate(density = V7, type1 = 'gene',
           type2 = 'protein_coding gene') %>%
    dplyr::select(chr = V1, bin = V4,
                  density, type1, type2),
  read.table(
    file = "34.LTR_chr//Aly/chr_sw.Copia.density.txt") %>%
    mutate(density = V7, type1 = 'TE',
           type2 = 'LTR/Copia') %>%
    dplyr::select(chr = V1, bin = V4,
                  density, type1, type2),
  read.table(
    file = "34.LTR_chr//Aly/chr_sw.Gypsy.density.txt") %>%
    mutate(density = V7, type1 = 'TE',
           type2 = 'LTR/Gypsy') %>%
    dplyr::select(chr = V1, bin = V4,
                  density, type1, type2)
) %>%
  mutate(Species = "Aly")
# filter(chr %in% 1:5) # in Ath

# distribution along the chromosome
Aly.chr_sw_density_summ = group_by(Aly.chr_sw_density, chr) %>%
  summarise(
    start = min(bin),
    end = max(bin)) %>%
  mutate(middle = start + (end - start)*0.5)
library(ggsci)
ggplot(data = Aly.chr_sw_dencity, aes(x = bin, y = density)) +
  geom_line(aes(color = type2), size = 0.1) +
  geom_vline(xintercept = Aly.chr_sw_density_summ$end,
             linetype = "dotdash") +
  facet_wrap(~type1,
             ncol = 1,
             scales = "free") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = Aly.chr_sw_density_summ$middle,
                     labels = Aly.chr_sw_density_summ$chr) +
  scale_color_lancet() +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_classic() +
    theme(
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0, 'lines'),
      strip.text = element_blank(),
      legend.position = "top"
    )

# same for Ath