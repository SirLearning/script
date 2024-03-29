library(tidyverse)
library(rtracklayer)

Aly.genome_size = sum(read.table(file = "../data/chr1A.fa.fai")$V2)
Ath.genome_size = sum(read.table(file = "../data/chr1A.fa.fai")$V2)

Aly.allTE = import("../data/chr1A.anno.gff3") %>%
  as_tibble() %>%
  mutate(Species = "Aly", genome_size = Aly.genome_size)

Ath.allTE = import("../data/filter.anno.gff3") %>%
  as_tibble() %>%
  mutate(Species = "Ath", genome_size = Ath.genome_size)

allTE = rbind(Aly.allTE, Ath.allTE) %>%
  mutate(Classification0 =
         str_split(Classification, "/", simplify = TRUE)[,1])
rm(Aly.allTE, Ath.allTE)

allTE_summ = filter(allTE,
                    !(type %in% c("repeat_region", "target_site_duplication", "long_terminal_repeat"))) %>%
  group_by(Species, Classification) %>%
    summarise(
      count = n(),
      mean_size = mean(width),
      size = sum(width),
      percent = sum(width) / mean(genome_size) * 100
    )
# stats of number of TEs
intactTE = filter(allTE, Method=="structural")
intactTE_summ = group_by(allTE, Species, Classification) %>%
  summarise(
    count = n(),
    size = sum(width),
    percent = sum(width) / mean(genome_size) * 100
  )
