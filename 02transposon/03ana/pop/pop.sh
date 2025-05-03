cat *lib.fa | grep ">" | grep -v "Unknown" | grep -v "unknown" > list
cut -c 2- list > list.cut
rm list && mv list.cut list
seqkit grep -f list ab.durum.fa.mod.EDTA.TElib.fa > ab.durum.known.fa
# check same seqs in libs
RepeatMasker -e ncbi -pa 10 -q -no_is -norna -nolow -div 20 -lib $TE2 $TE1 # 80% same
# we
sh ~/script/04ana/pop/04sample.sh CRR072247 ab.we
sh ~/script/04ana/pop/04sample.sh CRR072248 ab.we
sh ~/script/04ana/pop/04sample.sh CRR072249 ab.we
sh ~/script/04ana/pop/04sample.sh CRR072250 ab.we
sh ~/script/04ana/pop/04sample.sh CRR072251 ab.we
# dw
sh ~/script/04ana/pop/04sample.sh CRR072337 ab.durum
sh ~/script/04ana/pop/04sample.sh CRR072338 ab.durum
sh ~/script/04ana/pop/04sample.sh CRR072340 ab.durum
sh ~/script/04ana/pop/04sample.sh CRR072341 ab.durum
sh ~/script/04ana/pop/04sample.sh CRR072347 ab.durum
# at
sh ~/script/04ana/pop/04sample.sh CRR072405 d.tauchii
sh ~/script/04ana/pop/04sample.sh CRR072406 d.tauchii
sh ~/script/04ana/pop/04sample.sh CRR072407 d.tauchii
sh ~/script/04ana/pop/04sample.sh CRR072408 d.tauchii
sh ~/script/04ana/pop/04sample.sh CRR072409 d.tauchii
# ldr
sh ~/script/04ana/pop/03map_threshold.sh CRR072401 1A
sh ~/script/04ana/pop/03map_threshold.sh SRR7164572 1A
sh ~/script/04ana/pop/03map_threshold.sh SRR7164576 1A
sh ~/script/04ana/pop/03map_threshold.sh SRR7164580 1A
sh ~/script/04ana/pop/03map_threshold.sh SRR7164606 1A

sh ~/script/04ana/pop/03map_threshold.sh CRR072401 1B
sh ~/script/04ana/pop/03map_threshold.sh SRR7164572 1B
sh ~/script/04ana/pop/03map_threshold.sh SRR7164576 1B
sh ~/script/04ana/pop/03map_threshold.sh SRR7164580 1B
sh ~/script/04ana/pop/03map_threshold.sh SRR7164606 1B

sh ~/script/04ana/pop/03map_threshold.sh CRR072401 1D
sh ~/script/04ana/pop/03map_threshold.sh SRR7164572 1D
sh ~/script/04ana/pop/03map_threshold.sh SRR7164576 1D
sh ~/script/04ana/pop/03map_threshold.sh SRR7164580 1D
sh ~/script/04ana/pop/03map_threshold.sh SRR7164606 1D
# ctv
sh ~/script/04ana/pop/03map_threshold.sh SRR7164604 a.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164620 a.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164628 a.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164669 a.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164670 a.lancer

sh ~/script/04ana/pop/03map_threshold.sh SRR7164604 b.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164620 b.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164628 b.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164669 b.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164670 b.lancer

sh ~/script/04ana/pop/03map_threshold.sh SRR7164604 d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164620 d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164628 d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164669 d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164670 d.lancer