library(tidyverse)
library(rtracklayer)

Aly.genome_size = sum(read.table(file = "12.ref/Aly.genome.fa.fai")$V2)
Ath.genome_size = sum(read.table(file = "12.ref/Ath.genome.fa.fai")$V2)

Aly.allTE = import("23.EDTA/Aly/genome.fa.mod.EDTA.TEanno.gff3") %>%
  as_tibble() %>%
  mutate(Species = "Aly", genome_size = 1.5e8)

Ath.allTE = import("23.EDTA/Ath/genome.fa.mod.EDTA.TEanno.gff3") %>%
  as_tibble() %>%
  mutate(Species = "Ath", genome_size = 1.5e8)

allTE = rbind(Aly.allTE, Ath.allTE) %>%
  mutate(Classification0 =
         str_split(Classification, "/", simplify = TRUE)[,1])
rm(Aly.allTE, Ath.allTE)

allTE_summ = filter(allTE,
                    !(type %in% c("repeat_region", "target_site_duplication", "long_terminal_repeat"))) %>%
  group_by(Species, Classification0) %>%
    summarise(
      count = n(),
      mean_size = mean(width),
      size = sum(width),
      percent = sum(width) / mean(genome_size) * 100
    )