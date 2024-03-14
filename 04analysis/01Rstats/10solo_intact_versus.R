library(rtracklayer)
# extract all aligned LTR|INT annotations
LTR_RM = import(
  "23.EDTA/Aly/genome.fa.mod.EDTA.anno/genome.fa.mod.EDTA.RM.gff3") %>%
  as_tibble() %>%
  filter(str_detect(Name, "INT|LTR"))

library(readr)
# read the sequence length in LTRlib
LTRlib_size = read_table(
  "LTR_solo_intact_ratio/Aly/genome.fa.mod.EDTA.TElib.fa.fai",
    col_names = FALSE) %>%
  filter(str_detect(X1, "LTR|INT")) %>%
  separate(X1, sep = "#", into = c("Name", "Classification")) %>%
  select(Name, size = X2)

# extract the location of the LTR and INT separately
#   filter the alignment result
LTRRT_homo = inner_join(LTR_RM, LTRlib_size, by = "Name") %>%
  # score > 300, coverage > 80%
  filter(score > 300, width/size > 0.8) %>%
  mutate(region = str_sub(Name, 13, 15)) %>%
    select(sequence, start, end, ID, Name, region)
#   LTR alignment
LTR_homo = filter(LTRRT_homo, region == "LTR") %>%
  mutate(row_number = row_number())
#   INT alignment
INT_homo = filter(LTRRT_homo, region == "INT") %>%
  mutate(row_number = row_number())
#   LTR that is nearest to INT
library(hiAnnotator)
INT_Nearest_LTR = get2NearestFeature(
  sites.rd = as(INT_homo, "GRanges"),
  features.rd = as(LTR_homo, "GRanges"),
  "NearestLTR",
  side = "either") %>%
  as_tibble()

# extract intactLTR and soloLTR
#   intact LTRRT
ltr_int_dist = 1000 # both contain a LTR in the 100bp of two sides
intactLTR_RM =
  filter(INT_Nearest_LTR,
         Either.NearestLTR.upStream1.Dist >= -1*ltr_int_dist &
         Either.NearestLTR.downStrean1.Disr <= ltr_int_dist &
         Either.NearestLTR.upStream1 == Either.NearestLTR.downStream1)
#   solo LTRRT
soloLTR_RM =
  filter(INT_Nearest_LTR,
         (Either.NearestLTR.upStream1.Dist >= -1*ltr_int_dist |
           Either.NearestLTR.downStream1.Dist <= ltr_int_dist) &
           !(ID %in% intactLTR_RM$ID))