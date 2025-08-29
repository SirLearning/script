# we
nohup sh ~/script/WeaTE/pop/abd.sh CRR072247 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072248 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072249 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072250 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072251 we ab &

# dw
nohup sh ~/script/WeaTE/pop/abd.sh CRR072337 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072338 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072340 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072341 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072347 dw ab &

# at
nohup sh ~/script/WeaTE/pop/abd.sh CRR072405 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072406 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072407 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072408 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072409 at d &


# cs
sh ~/script/WeaTE/pop/abd.sh CRR072401 cs a
sh ~/script/WeaTE/pop/abd.sh CRR072401 cs b
cat a.CRR072401.sum b.CRR072401.sum > ab.CRR072401.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164576 cs a
sh ~/script/WeaTE/pop/abd.sh SRR7164576 cs b
cat a.SRR7164576.sum b.SRR7164576.sum > ab.SRR7164576.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164580 cs a
sh ~/script/WeaTE/pop/abd.sh SRR7164580 cs b
cat a.SRR7164580.sum b.SRR7164580.sum > ab.SRR7164580.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164572 cs a
sh ~/script/WeaTE/pop/abd.sh SRR7164572 cs b
cat a.SRR7164572.sum b.SRR7164572.sum > ab.SRR7164572.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164606 cs a
sh ~/script/WeaTE/pop/abd.sh SRR7164606 cs b
cat a.SRR7164606.sum b.SRR7164606.sum > ab.SRR7164606.sum

nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 cs d &

# la
sh ~/script/WeaTE/pop/abd.sh SRR7164620 la a
sh ~/script/WeaTE/pop/abd.sh SRR7164620 la b
cat a.SRR7164620.sum b.SRR7164620.sum > ab.SRR7164620.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164628 la a
sh ~/script/WeaTE/pop/abd.sh SRR7164628 la b
cat a.SRR7164628.sum b.SRR7164628.sum > ab.SRR7164628.sumd ~/script/WeaTE/pop/abd.sh SRR7164670 la a
sh ~/script/WeaTE/pop/abd.sh SRR7164670 la a
sh ~/script/WeaTE/pop/abd.sh SRR7164670 la b
cat a.SRR7164670.sum b.SRR7164670.sum > ab.SRR7164670.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164669 la a
sh ~/script/WeaTE/pop/abd.sh SRR7164669 la b
cat a.SRR7164669.sum b.SRR7164669.sum > ab.SRR7164669.sum
sh ~/script/WeaTE/pop/abd.sh SRR7164604 la a
sh ~/script/WeaTE/pop/abd.sh SRR7164604 la b
cat a.SRR7164604.sum b.SRR7164604.sum > ab.SRR7164604.sum

nohup sh ~/script/WeaTE/pop/abd.sh SRR7164620 la d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164628 la d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164670 la d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164669 la d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164604 la d &



# !!!! CS !!!!
# 1A sh
sh ~/script/WeaTE/pop/cs.sh a bcs cs
sh ~/script/WeaTE/pop/cs.sh a dcs cs
sh ~/script/WeaTE/pop/cs.sh a dw cs
sh ~/script/WeaTE/pop/cs.sh a we cs
mv a-*cs*.sum ..

# 1B sh
sh ~/script/WeaTE/pop/cs.sh b dw cs
sh ~/script/WeaTE/pop/cs.sh b we cs
nohup sh ~/script/WeaTE/pop/cs.sh b acs cs &
nohup sh ~/script/WeaTE/pop/cs.sh b dcs cs &
mv b-*cs*.sum ..

# 1D sh
sh ~/script/WeaTE/pop/cs.sh d at cs
nohup sh ~/script/WeaTE/pop/cs.sh d acs cs &
nohup sh ~/script/WeaTE/pop/cs.sh d bcs cs &

# .. sh
nohup cat a-bw.CRR072401.sum b-bw.CRR072401.sum > ab-bw.CRR072401.sum
nohup cat a-bw.SRR7164576.sum b-bw.SRR7164576.sum > ab-bw.SRR7164576.sum
nohup cat a-bw.SRR7164580.sum b-bw.SRR7164580.sum > ab-bw.SRR7164580.sum
nohup cat a-bw.SRR7164572.sum b-bw.SRR7164572.sum > ab-bw.SRR7164572.sum
nohup cat a-bw.SRR7164606.sum b-bw.SRR7164606.sum > ab-bw.SRR7164606.sum

cat a-we.CRR072401.sum b-we.CRR072401.sum > ab-we.CRR072401.sum
cat a-we.SRR7164576.sum b-we.SRR7164576.sum > ab-we.SRR7164576.sum
cat a-we.SRR7164580.sum b-we.SRR7164580.sum > ab-we.SRR7164580.sum
cat a-we.SRR7164572.sum b-we.SRR7164572.sum > ab-we.SRR7164572.sum
cat a-we.SRR7164606.sum b-we.SRR7164606.sum > ab-we.SRR7164606.sum



# !!!! new TE !!!!
nohup sh ~/script/WeaTE/pop/cs.sh a 1 ol &
nohup sh ~/script/WeaTE/pop/cs.sh a 2 ol &
nohup sh ~/script/WeaTE/pop/cs.sh a 3 ol &

nohup sh ~/script/WeaTE/pop/cs.sh b 1 ol &
nohup sh ~/script/WeaTE/pop/cs.sh b 2 ol &
nohup sh ~/script/WeaTE/pop/cs.sh b 3 ol &

nohup sh ~/script/WeaTE/pop/cs.sh d 1 ol &
nohup sh ~/script/WeaTE/pop/cs.sh d 2 ol &

nohup sh ~/script/WeaTE/pop/cs.sh a 0 ol &



# !!!! WE domestication !!!!
sh ~/script/WeaTE/pop/rm.sh we at # we-at
sh ~/script/WeaTE/pop/rm.sh we dw # we-dw
sh ~/script/WeaTE/pop/rm.sh dw at # dw-at
sh ~/script/WeaTE/pop/rm.sh dw we # dw-we

cat /data1/home/dazheng/transposon/pop/03mapping/ab.we | grep -v -f we-at > we-0 # ancestor-WE
cat we-at | grep -v -f we-dw > we-1 # WE-DW
cat we-at | grep -f we-dw > we-2 # WE-now

cat /data1/home/dazheng/transposon/pop/03mapping/ab.dw | grep -v -f dw-at > dw-0 # ancestor-DW
cat dw-at | grep -v -f dw-we > dw-1 # DW-WE
cat dw-at | grep -f dw-we > dw-2 # DW-now


# we
nohup sh ~/script/WeaTE/pop/abd.sh CRR072247 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072248 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072249 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072250 we ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072251 we ab &
# dw
nohup sh ~/script/WeaTE/pop/abd.sh CRR072337 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072338 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072340 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072341 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072347 dw ab &


# !!!!!! polyploidy !!!!!!
# a
cat /data1/home/dazheng/transposon/pop/03mapping/a.cs | grep -v -f cs-we > cs-0 # ancestor-DW
cat cs-we | grep -v -f cs-dw > cs-1 # DW-CS
cat cs-we | grep -f cs-dw > cs-2 # CS-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 cs a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 cs a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 cs a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 cs a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 cs a &
# b
cat /data1/home/dazheng/transposon/pop/03mapping/b.cs | grep -v -f cs-we > cs-0 # ancestor-DW
cat cs-we | grep -v -f cs-dw > cs-1 # DW-CS
cat cs-we | grep -f cs-dw > cs-2 # CS-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 cs b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 cs b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 cs b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 cs b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 cs b &

cat 01A/cs-0.CRR072401.sum 02B/cs-0.CRR072401.sum > 04AB/cs-0.CRR072401.sum
cat 01A/cs-0.SRR7164576.sum 02B/cs-0.SRR7164576.sum > 04AB/cs-0.SRR7164576.sum
cat 01A/cs-0.SRR7164580.sum 02B/cs-0.SRR7164580.sum > 04AB/cs-0.SRR7164580.sum
cat 01A/cs-0.SRR7164572.sum 02B/cs-0.SRR7164572.sum > 04AB/cs-0.SRR7164572.sum
cat 01A/cs-0.SRR7164606.sum 02B/cs-0.SRR7164606.sum > 04AB/cs-0.SRR7164606.sum
cat 01A/cs-1.CRR072401.sum 02B/cs-1.CRR072401.sum > 04AB/cs-1.CRR072401.sum
cat 01A/cs-1.SRR7164576.sum 02B/cs-1.SRR7164576.sum > 04AB/cs-1.SRR7164576.sum
cat 01A/cs-1.SRR7164580.sum 02B/cs-1.SRR7164580.sum > 04AB/cs-1.SRR7164580.sum
cat 01A/cs-1.SRR7164572.sum 02B/cs-1.SRR7164572.sum > 04AB/cs-1.SRR7164572.sum
cat 01A/cs-1.SRR7164606.sum 02B/cs-1.SRR7164606.sum > 04AB/cs-1.SRR7164606.sum
cat 01A/cs-2.CRR072401.sum 02B/cs-2.CRR072401.sum > 04AB/cs-2.CRR072401.sum
cat 01A/cs-2.SRR7164576.sum 02B/cs-2.SRR7164576.sum > 04AB/cs-2.SRR7164576.sum
cat 01A/cs-2.SRR7164580.sum 02B/cs-2.SRR7164580.sum > 04AB/cs-2.SRR7164580.sum
cat 01A/cs-2.SRR7164572.sum 02B/cs-2.SRR7164572.sum > 04AB/cs-2.SRR7164572.sum
cat 01A/cs-2.SRR7164606.sum 02B/cs-2.SRR7164606.sum > 04AB/cs-2.SRR7164606.sum


# dw
sh ~/script/WeaTE/pop/rm.sh dw acs # dw-acs
sh ~/script/WeaTE/pop/rm.sh dw bcs # dw-bcs
cat dw-acs | grep -f dw-bcs > dw-cs
sh ~/script/WeaTE/pop/rm.sh dw we

cat /data1/home/dazheng/transposon/pop/03mapping/ab.dw | grep -v -f dw-we > dw-0 # ancestor-DW
cat dw-we | grep -v -f dw-cs > dw-1 # DW-CS
cat dw-we | grep -f dw-cs > dw-2 # DW-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072337 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072338 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072340 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072341 dw ab &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072347 dw ab &


# d
cat /data1/home/dazheng/transposon/pop/03mapping/d.cs | grep -v -f cs-acs > cs-0 # ancestor-AT
cat cs-acs | grep -v -f cs-at > cs-1 # AT-CS
cat cs-acs | grep -f cs-at > cs-2 # CS-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 cs d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 cs d &


# at
sh ~/script/WeaTE/pop/rm.sh at acs # at-acs
sh ~/script/WeaTE/pop/rm.sh at cs # at-cs

cat /data1/home/dazheng/transposon/pop/03mapping/d.at | grep -v -f at-acs > at-0 # ancestor-AT
cat at-acs | grep -v -f at-cs > at-1 # AT-CS
cat at-acs | grep -f at-cs > at-2 # AT-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072405 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072406 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072407 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072408 at d &
nohup sh ~/script/WeaTE/pop/abd.sh CRR072409 at d &


# !!!!!! BW domestication !!!!!!
# a
sh ~/script/WeaTE/pop/rm.sh la dw # la-dw
sh ~/script/WeaTE/pop/rm.sh la cs # la-cs
cat /data1/home/dazheng/transposon/pop/03mapping/a.cs | grep -v -f cs-dw > ldr-0 # ancestor-landrace
cat cs-dw | grep -v -f cs-la > ldr-1 # landrace-cultivar
cat cs-dw | grep -f cs-la > ldr-2 # landrace-now
cat /data1/home/dazheng/transposon/pop/03mapping/a.la | grep -v -f la-dw > ctv-0 # ancestor-landrace
cat la-dw | grep -v -f la-cs > ctv-1 # landrace-cultivar
cat la-dw | grep -f la-cs > ctv-2 # cultivar-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 ldr a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 ldr a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 ldr a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 ldr a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 ldr a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164620 ctv a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164628 ctv a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164670 ctv a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164669 ctv a &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164604 ctv a &
# b
sh ~/script/WeaTE/pop/rm.sh la dw # la-dw
sh ~/script/WeaTE/pop/rm.sh la cs # la-cs
cat /data1/home/dazheng/transposon/pop/03mapping/b.cs | grep -v -f cs-dw > ldr-0 # ancestor-landrace
cat cs-dw | grep -v -f cs-la > ldr-1 # landrace-cultivar
cat cs-dw | grep -f cs-la > ldr-2 # landrace-now
cat /data1/home/dazheng/transposon/pop/03mapping/b.la | grep -v -f la-dw > ctv-0 # ancestor-landrace
cat la-dw | grep -v -f la-cs > ctv-1 # landrace-cultivar
cat la-dw | grep -f la-cs > ctv-2 # cultivar-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 ldr b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 ldr b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 ldr b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 ldr b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 ldr b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164620 ctv b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164628 ctv b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164670 ctv b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164669 ctv b &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164604 ctv b &
# d
sh ~/script/WeaTE/pop/rm.sh la at # la-at
sh ~/script/WeaTE/pop/rm.sh la cs # la-cs
cat /data1/home/dazheng/transposon/pop/03mapping/d.cs | grep -v -f cs-at > ldr-0 # ancestor-landrace
cat cs-at | grep -v -f cs-la > ldr-1 # landrace-cultivar
cat cs-at | grep -f cs-la > ldr-2 # landrace-now
cat /data1/home/dazheng/transposon/pop/03mapping/d.la | grep -v -f la-at > ctv-0 # ancestor-landrace
cat la-at | grep -v -f la-cs > ctv-1 # landrace-cultivar
cat la-at | grep -f la-cs > ctv-2 # cultivar-now
nohup sh ~/script/WeaTE/pop/abd.sh CRR072401 ldr d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164576 ldr d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164580 ldr d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164572 ldr d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164606 ldr d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164620 ctv d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164628 ctv d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164670 ctv d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164669 ctv d &
nohup sh ~/script/WeaTE/pop/abd.sh SRR7164604 ctv d &


cat 01A/ldr-0.CRR072401.sum 02B/ldr-0.CRR072401.sum 03D/ldr-0.CRR072401.sum > ldr-0.CRR072401.sum
cat 01A/ldr-0.SRR7164576.sum 02B/ldr-0.SRR7164576.sum 03D/ldr-0.SRR7164576.sum > ldr-0.SRR7164576.sum
cat 01A/ldr-0.SRR7164580.sum 02B/ldr-0.SRR7164580.sum 03D/ldr-0.SRR7164580.sum > ldr-0.SRR7164580.sum
cat 01A/ldr-0.SRR7164572.sum 02B/ldr-0.SRR7164572.sum 03D/ldr-0.SRR7164572.sum > ldr-0.SRR7164572.sum
cat 01A/ldr-0.SRR7164606.sum 02B/ldr-0.SRR7164606.sum 03D/ldr-0.SRR7164606.sum > ldr-0.SRR7164606.sum

cat 01A/ldr-1.CRR072401.sum 02B/ldr-1.CRR072401.sum 03D/ldr-1.CRR072401.sum > ldr-1.CRR072401.sum
cat 01A/ldr-1.SRR7164576.sum 02B/ldr-1.SRR7164576.sum 03D/ldr-1.SRR7164576.sum > ldr-1.SRR7164576.sum
cat 01A/ldr-1.SRR7164580.sum 02B/ldr-1.SRR7164580.sum 03D/ldr-1.SRR7164580.sum > ldr-1.SRR7164580.sum
cat 01A/ldr-1.SRR7164572.sum 02B/ldr-1.SRR7164572.sum 03D/ldr-1.SRR7164572.sum > ldr-1.SRR7164572.sum
cat 01A/ldr-1.SRR7164606.sum 02B/ldr-1.SRR7164606.sum 03D/ldr-1.SRR7164606.sum > ldr-1.SRR7164606.sum

cat 01A/ldr-2.CRR072401.sum 02B/ldr-2.CRR072401.sum 03D/ldr-2.CRR072401.sum > ldr-2.CRR072401.sum
cat 01A/ldr-2.SRR7164576.sum 02B/ldr-2.SRR7164576.sum 03D/ldr-2.SRR7164576.sum > ldr-2.SRR7164576.sum
cat 01A/ldr-2.SRR7164580.sum 02B/ldr-2.SRR7164580.sum 03D/ldr-2.SRR7164580.sum > ldr-2.SRR7164580.sum
cat 01A/ldr-2.SRR7164572.sum 02B/ldr-2.SRR7164572.sum 03D/ldr-2.SRR7164572.sum > ldr-2.SRR7164572.sum
cat 01A/ldr-2.SRR7164606.sum 02B/ldr-2.SRR7164606.sum 03D/ldr-2.SRR7164606.sum > ldr-2.SRR7164606.sum

cat 01A/ctv-0.SRR7164620.sum 02B/ctv-0.SRR7164620.sum 03D/ctv-0.SRR7164620.sum > ctv-0.SRR7164620.sum
cat 01A/ctv-0.SRR7164628.sum 02B/ctv-0.SRR7164628.sum 03D/ctv-0.SRR7164628.sum > ctv-0.SRR7164628.sum
cat 01A/ctv-0.SRR7164670.sum 02B/ctv-0.SRR7164670.sum 03D/ctv-0.SRR7164670.sum > ctv-0.SRR7164670.sum
cat 01A/ctv-0.SRR7164669.sum 02B/ctv-0.SRR7164669.sum 03D/ctv-0.SRR7164669.sum > ctv-0.SRR7164669.sum
cat 01A/ctv-0.SRR7164604.sum 02B/ctv-0.SRR7164604.sum 03D/ctv-0.SRR7164604.sum > ctv-0.SRR7164604.sum

cat 01A/ctv-1.SRR7164620.sum 02B/ctv-1.SRR7164620.sum 03D/ctv-1.SRR7164620.sum > ctv-1.SRR7164620.sum
cat 01A/ctv-1.SRR7164628.sum 02B/ctv-1.SRR7164628.sum 03D/ctv-1.SRR7164628.sum > ctv-1.SRR7164628.sum
cat 01A/ctv-1.SRR7164670.sum 02B/ctv-1.SRR7164670.sum 03D/ctv-1.SRR7164670.sum > ctv-1.SRR7164670.sum
cat 01A/ctv-1.SRR7164669.sum 02B/ctv-1.SRR7164669.sum 03D/ctv-1.SRR7164669.sum > ctv-1.SRR7164669.sum
cat 01A/ctv-1.SRR7164604.sum 02B/ctv-1.SRR7164604.sum 03D/ctv-1.SRR7164604.sum > ctv-1.SRR7164604.sum

cat 01A/ctv-2.SRR7164620.sum 02B/ctv-2.SRR7164620.sum 03D/ctv-2.SRR7164620.sum > ctv-2.SRR7164620.sum
cat 01A/ctv-2.SRR7164628.sum 02B/ctv-2.SRR7164628.sum 03D/ctv-2.SRR7164628.sum > ctv-2.SRR7164628.sum
cat 01A/ctv-2.SRR7164670.sum 02B/ctv-2.SRR7164670.sum 03D/ctv-2.SRR7164670.sum > ctv-2.SRR7164670.sum
cat 01A/ctv-2.SRR7164669.sum 02B/ctv-2.SRR7164669.sum 03D/ctv-2.SRR7164669.sum > ctv-2.SRR7164669.sum
cat 01A/ctv-2.SRR7164604.sum 02B/ctv-2.SRR7164604.sum 03D/ctv-2.SRR7164604.sum > ctv-2.SRR7164604.sum

