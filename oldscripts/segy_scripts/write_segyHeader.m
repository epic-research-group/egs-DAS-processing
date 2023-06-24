clear;clc;
!ls /data1/parker/EGS_CASSM2/* -d > ff.dat
!wc -l ff.dat > nd.txt
fid=fopen('nd.txt');
nd=textscan(fid,'%d %s');
nd=nd{1};
i=1;
f=regexp(fileread('ff.dat'),'\n','split');
for i=1:nd
    [fol,v1,~]=fileparts(f{i});
    cd([fol,'/',v1])
    if isempty(dir('*.dat'))
        continue   
    end
    if ~isempty(dir('20*.log'))
        !grep -c ".dat" 20*.log | grep -n ".dat" 20*.log | cut -d" " -f3 > ns.txt
        ns=load('ns.txt'); ns=ns+1;
        
        gf='!cat ';
        for j=1:length(ns)
            %gf=cat(1,gf,load(['~/scripts/segy_scripts/',num2str(ns(j)),'.dat']));
            gf=[gf ['~/scripts/segy_scripts/',num2str(ns(j)),'.dat ']];
            %should I write this as a .dat file in this directory?
        end
        gf=[gf '> gf.dat'];
        eval(gf);
        try
            !rm geom*.sgy
        catch me
        end
        
        sg=dir('*.sgy');
        sg=sg(1).name;
        eval(['!segyread verbose=0 tape=' sg ' > out.su'])
        !a2b verbose=0 < gf.dat > geom.bin
        
        !sushw verbose=0 < out.su infile=geom.bin key=sx,sy,selev,gx,gy,gelev,scalel,scalco,dt > out1.su
        
        !segyhdrs verbose=0 < out1.su
        
        eval(['!segywrite verbose=0 tape=geom_' sg ' < out1.su'])
        
        !rm geom.bin
        !rm  gf.dat 
        !rm *.su
    else
        continue
    end
    i = i+1;
    
end
