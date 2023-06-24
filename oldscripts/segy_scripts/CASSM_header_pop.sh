#! /bin/sh
# populate headers in the SEGY files for the EGS Collab CASSM data
ls /data1/parker/EGS_CASSM/* -d > ff.dat

wc -l ff.dat > nd.txt
nd=`awk '{print $1}' nd.txt`
i=1
while [ $i -le $nd ]
do
	v1=`awk '(NR=='$i'){printf $1}' /home/spri902/scripts/segy_scripts/ff.dat`
	cd $v1
	nf=grep -c ".dat" 20*.log
	if [ $nf -lt 17 ]

	then

	grep -c ".dat" 20*V.log | grep -n ".dat" 20*.log | cut -d" " -f3 > ns.txt
	cat 	
	else
		
	ls *sgy > s.txt
	f1=`awk '(NR==1){printf $1}' s.txt
	/usr/local/cwp/src/Third_Party/segyread/segyread tape=$f1 > out.su
	
	a2b < geomfile.dat > geom.bin

	sushw < out.su infile=geom.bin key=sx,sy,selev,gx,gy,gelev,scalel,scalco,dt > out1.su

	segyhdrs < out1.su

	segywrite tape=geom_$f1.sgy < out1.su

	rm stuff
























