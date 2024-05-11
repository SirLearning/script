# a
cat cs-bcs | grep -f cs-dcs > no-bd # 2674
cat /data1/home/dazheng/transposon/pop/03mapping/a.cs | grep -v -f no-bd > ol-0 # 624
cat no-bd | grep -f cs-dw | grep -f cs-we > ol-cs # 1456
cat no-bd | grep -v -f ol-cs > no-bd-cs # 1218
cat no-bd-cs | grep -f cs-we > ol-dw # 317
cat no-bd-cs | grep -v -f ol-dw > ol-we # 901




# b
cat cs-acs | grep -f cs-dcs > no-ad
cat /data1/home/dazheng/transposon/pop/03mapping/b.cs | grep -v -f no-ad > ol-0
cat no-ad | grep -f cs-dw | grep -f cs-we > ol-cs
cat no-ad | grep -v -f ol-cs > no-ad-cs
cat no-ad-cs | grep -f cs-we > ol-dw
cat no-ad-cs | grep -v -f ol-dw > ol-we

mv ol-cs ol-3
mv ol-dw ol-2
mv ol-we ol-1

# d
cat cs-acs | grep -f cs-bcs > no-ab
cat /data1/home/dazheng/transposon/pop/03mapping/d.cs | grep -v -f no-ab > ol-0
cat no-ab | grep -f cs-at > ol-cs
cat no-ab | grep -v -f ol-cs > ol-at

mv ol-cs ol-2
mv ol-at ol-1