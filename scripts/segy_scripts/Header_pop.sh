#!/bin/sh


segyread tape=filename.sgy > out.su

awk '{printf columns that you care about}' geomfile > geominput.dat
a2b < geominput.dat > geom.bin

sushw < out.su infile=geom.bin key=sx,sy,selev,rx,ry,gelev,scacleco,scalel > out1.su

segyhdrs < out1.su

segywrite tape=segyfilewithgeom.sgy < out1.su


