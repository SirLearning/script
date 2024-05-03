# we
sh ~/script/04ana/pop/03map_threshold.sh CRR072337 ab.durum
sh ~/script/04ana/pop/03map_threshold.sh CRR072338 ab.durum
sh ~/script/04ana/pop/03map_threshold.sh CRR072340 ab.durum
sh ~/script/04ana/pop/03map_threshold.sh CRR072341 ab.durum
sh ~/script/04ana/pop/03map_threshold.sh CRR072347 ab.durum

#sh ~/script/04ana/pop/03map_threshold.sh CRR072405 d.tauchii
#sh ~/script/04ana/pop/03map_threshold.sh CRR072406 d.tauchii
#sh ~/script/04ana/pop/03map_threshold.sh CRR072407 d.tauchii
#sh ~/script/04ana/pop/03map_threshold.sh CRR072408 d.tauchii
#sh ~/script/04ana/pop/03map_threshold.sh CRR072409 d.tauchii

cat *lib.fa | grep ">" | grep -v "Unknown" | grep -v "unknown" > list
cut -c 2- list > list.cut
rm list && mv list.cut list
seqkit grep -f list ab.durum.fa.mod.EDTA.TElib.fa > ab.durum.known.fa

#nohup sh ~/script/04ana/pop/03map_threshold.sh SRR7164580 1A &
#nohup sh ~/script/04ana/pop/03map_threshold.sh SRR7164606 1A &

sh ~/script/04ana/pop/03map_threshold.sh SRR7164604 01A a.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164604 02B b.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164604 03D d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164620 01A a.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164620 02B b.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164620 03D d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164628 01A a.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164628 02B b.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164628 03D d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164669 01A a.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164669 02B b.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164669 03D d.lancer
sh ~/script/04ana/pop/03map_threshold.sh SRR7164670 01A a.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164670 02B b.lancer ; sh ~/script/04ana/pop/03map_threshold.sh SRR7164670 03D d.lancer








RepeatMasker -e ncbi -pa 10 -q -no_is -norna -nolow -div 20 -lib $TE2 $TE1 # 80% same
