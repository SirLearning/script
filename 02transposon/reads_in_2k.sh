#awk '/^>/{if(NR>1){x++}if(x%2000==0){f=sprintf("file%d.fasta",x)};print > f;next}{print > f}' ../1b.2x.clean.fasta
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file0.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog0 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file2000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog2000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file4000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog4000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file6000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog6000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file8000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog8000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file10000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog10000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file12000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog12000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file14000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog14000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file16000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog16000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file18000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog18000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file20000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog20000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file22000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog22000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file24000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog24000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file26000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog26000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file28000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog28000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file30000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog30000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file32000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog32000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file34000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog34000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file36000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog36000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file38000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog38000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file40000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog40000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file42000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog42000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file44000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog44000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file46000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog46000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file48000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog48000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file50000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog50000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file52000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog52000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file54000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog54000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file56000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog56000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file58000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog58000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file60000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog60000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file62000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog62000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file64000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog64000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file66000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 20 > alog66000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file68000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog68000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file70000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog70000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file72000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog72000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file74000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog74000 2>&1 &
sleep 18000
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file76000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog76000 2>&1 &
nohup perl ~/transposon/tools/EDTA/EDTA.pl --genome file78000.fasta --curatedlib ~/transposon/ABD/00data/00library/00curatedLib/abd.fa --overwrite 0 --sensitive 1 --threads 10 > alog78000 2>&1 &