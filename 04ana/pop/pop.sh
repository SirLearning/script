# we
sh ~/script/04ana/pop/03map_threshold.sh CRR072337 ab.we
sh ~/script/04ana/pop/03map_threshold.sh CRR072338 ab.we
sh ~/script/04ana/pop/03map_threshold.sh CRR072340 ab.we
sh ~/script/04ana/pop/03map_threshold.sh CRR072341 ab.we
sh ~/script/04ana/pop/03map_threshold.sh CRR072347 ab.we

#cat *lib.fa | grep ">" | grep -v "Unknown" | grep -v "unknown" > list
#cut -c 2- list > list.cut
#rm list && mv list.cut list
#seqkit grep -f list ab.durum.fa.mod.EDTA.TElib.fa > ab.durum.known.fa