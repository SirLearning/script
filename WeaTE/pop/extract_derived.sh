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
sh ~/script/WeaTE/pop/rm.sh dw we # dw-dw

cat /data1/home/dazheng/transposon/pop/03mapping/ab.we | grep -v -f we-at > we-0 # ancestor-WE
cat we-at | grep -v -f we-dw > we-1 # WE-DW
cat we-at | grep -f we-dw > we-2 # WE-now

cat /data1/home/dazheng/transposon/pop/03mapping/ab.dw | grep -v -f dw-at > dw-0 # ancestor-DW
cat dw-at | grep -v -f dw-dw > dw-1 # DW-WE
cat dw-at | grep -f dw-dw > dw-2 # DW-now


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

