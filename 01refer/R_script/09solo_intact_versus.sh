ln -s ../../23.EDTA/Aly/genome.fa.mod.EDTA.anno/genome.fa.mod.EDTA.RM.gff3 .
ln -s ../../23.EDTA/Aly/genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa .

# TElib index
samtools faidx genome.fa.mod.EDTA.TElib.fa

# extract solo & intact LTR
Rscript ../../11.software/solo_intact_ltr_edta.R \
  -g genome.fa.mod.EDTA.RM.gff3 # gff after RepeatMasker \
  -T genome.fa.mod.EDTA.TElib.fa.fai # using TElib as library