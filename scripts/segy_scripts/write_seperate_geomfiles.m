clear all; close all; clc;

src=load('srcfile.csv');
rec=load('recfile.csv');

sc=1000;
src=round(src*sc);
rec=round(rec*sc);
dt=48000;
for i=1:length(src(:,1))
    
    fid=fopen(sprintf('%i.dat',i),'w+');
    for j=1:length(rec(:,1))
    
    fprintf(fid,'%10d %10d %10d %10d %10d %10d %5d %5d %5d\n',...
        src(i,1),src(i,2),src(i,3),rec(j,1),rec(j,2),rec(j,3),sc,sc,dt);
    end
    fclose(fid);
end

