#! /bin/sh

# convert seg2 files to segy with SU Third Party seg2segy

ls /data1/parker/EGS_CASSM/* -d > ff2.dat

wc -l ff2.dat > nd2.txt
nd=`awk '{print $1}' nd2.txt`
i=1
while [ $i -le $nd ]

do 
echo $i
v1=`awk '(NR=='$i'){printf $1}' /data1/parker/EGS_CASSM/ff2.dat` 
echo $v1
cd $v1 
ls *dat > s2.txt 

f1=`awk '(NR==1){printf $1}' s2.txt` 

wc -l s2.txt > nf.txt

nf=`awk '{printf $1}'  nf.txt` 
#syntax is seg2segy first_file number_of_files
/usr/local/cwp/src/Third_Party/seg2segy/seg2segy $f1 $nf


	i=$(( $i + 1 ))

	rm s2.txt
	rm nf.txt
	rm tmp.txt

	
done
