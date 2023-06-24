addpath(genpath('/home/spri902/packages/segymat-1.6'))
addpath(genpath('/home/spri902/EGS_Collab/maystim'))
!ls 2*/geom_*.sgy >& list.txt
clear;clc;
%key=sx,sy,selev,gx,gy,gelev,scalel,scalco,dt


filename = 'list.txt';
A=importdata(filename);

for i=1:size(A,1)
    [Data{i},SegyTraceHeaders{i},SegyHeader{i}]=ReadSegy(A{i});
    
end
Data = Data(1:1000,:,:);
% srcfile=zeros(96,9,20);
% usrc = zeros(20,3);
% for i=1:20
%     srcfile(:,:,i)=load([num2str(i),'.dat']);
%     usrc(i,:) = srcfile(1,1:3,i);
% end  

[usx,sxi]=unique([SegyTraceHeaders.SourceX]);
[usy,syi]=unique([SegyTraceHeaders.SourceY]);
[uselev,selevi]=unique([SegyTraceHeaders.SourceSurfaceElevation]);
datFiles = unique([SegyTraceHeaders.FieldRecord]);
sortInd = sort(sxi);
for j = 1:size(Data,3)
    for i = 1:length(sortInd)-1
    
        src{i,j} = Data(:,sortInd(i):sortInd(i+1)-1,j);
  
    end
end    

for i=1:size(src,2)
    test12(:,i)=src{6,i}(:,12);
    test21(:,i)=src{6,i}(:,21);
    test22(:,i)=src{6,i}(:,22);
    test23(:,i)=src{6,i}(:,23);
    imagesc(src{6,i}(:,12:23))
end   
x22=max(max(abs(test22)));
for i=1:6
plot(test22(:,i)/x22+(i-1))
hold on
end