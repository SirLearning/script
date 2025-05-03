library(tidyverse)
library(rtracklayer)
# import the gff3 file
Aly.intactLTR_gff = import("31.intactLTR/Aly/intact.LTR.gff3") %>%
  as_tibble() %>%
  select(seqnames, start, end, width, strand, ID, Name, ltr_identity, motif, tsd)
# import insertion time
Aly.intactLTR_age = read_table(
  "23.EDTA/Aly/genome.fa.mod.EDTA.raw/LTR/genome.fa.mod.pass.list") %>%
  select(LTR_loc = `#LTR_loc`, Insertion_Time) %>%
  distinct()
# merge the two tables
Aly.intactLTR = mutate(
  Aly.intactLTR_gff,
  LTR_loc = str_c(seqnames, ":", start, "-", end)) %>%
  inner_join(Aly.intactLTR_age, by = "LTR_loc")
# read new classification
Aly.intactLTR_cls =
  read_delim("31.intactLTR/Aly/intact.LTR.fa.rexdb.cls.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# merge the classification
Aly.intactLTR = left_join(Aly.intactLTR, Aly.intactLTR_cls, by = c("ID" = "#TE")) %>%
  select(ID, Name, seqnames, start, end, width, strand, ltr_identity, motif, tsd, Order, Superfamily, Clade, Complete, Domains, Insertion_Time)

# same for Ath

# combine the two species
intactLTR = rbind(
  mutate(Aly.intactLTR, Species = "Aly"),
  mutate(Ath.intactLTR, Species = "Ath")
)
rm(Aly.intactLTR_age, Aly.intactLTR_cls, Aly.intactLTR_gff)

# insertion time estimation
filter(intactLTR, !is.na(Superfamily)) %>%
  ggplot(aes(x = Insertion_Time)) +
  geom_density(aes(y=..density..,
                   color=Species,
                   linetype=Superfamily)) +
  scale_x_continuous(
    name = "Insertion Time (Mys)",
    breaks = seq(0, 6000000, 500000), labels = seq(0, 6, 0.5)) +
  scale_color_lancet() +
  theme_classic() +
  theme(legend.position = c(0.8, 0.7))