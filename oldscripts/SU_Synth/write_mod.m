clear all; close all; clc;

dx=0.25; % x interval for model
dz=0.25; % z interval for model
m=1200;  % number of rows
n=1200;  % number of columns

sx=150;  % sourcex
gzx=180; % x location of vsp receivers

x=(0:n-1)*dx;
z=(0:m-1)*dz;
% read existing model
fid=fopen('c11_file','r');
c11=fread(fid,[m,n],'float32');
fclose(fid);
fid=fopen('c55_file','r');
c55=fread(fid,[m,n],'float32');
fclose(fid);
fid=fopen('rho_file','r');
rho=fread(fid,[m,n],'float32');
fclose(fid);

vp=sqrt(c11./rho);
vs=sqrt(c55./rho);



% or define model
% vp=4500*ones(m,n); % m/s
% vs=0.5*vp;
% rho=2500*ones(m,n); % kg/m^3
% vp(z<20,:)=2000;
% vs=0.5*vp;
% rho(z<20,:)=2000;

% add scatterers

I=find(z>51&z<=56);
for i=1:length(I)
    
    ind=randi(1200,1,300);
    vp(I(i),ind)=400;
    vs(I(i),ind)=0;
    
    rho(I(i),ind)=1600;
    
end

I=find(z>97&z<=100);

for i=1:length(I)
    
    ind=randi(1200,1,400);
    vp(I(i),ind)=1000;
    vs(I(i),ind)=500;
    
    rho(I(i),ind)=1600;
    
end



%% Plot model
figure;
subplot(121)
imagesc(x-150,z,vp);
hold on;
plot(0,41,'r*');
plot(0,55,'r*');
plot(0,61,'r*');
plot([30 30],[0 100],'k');
text(30,-5,'Well 9');
xlabel('meters');
ylabel('meters');
ylim([-10 100]);
xlim([-5; 40]);
ch=colorbar;
ylabel(ch,'Vp (m/s)');
set(gca,'fontsize',12,'fontweight','bold');

subplot(122)
imagesc(x-150,z,vs);
hold on;
plot(0,41,'r*');
plot(0,55,'r*');
plot(0,61,'r*');
plot([30 30],[0 100],'k');
text(30,-5,'Well 9');
xlabel('meters');
ylabel('meters');
ylim([-10 100]);
xlim([-5; 40]);
ch=colorbar;
ylabel(ch,'Vs (m/s)');
set(gca,'fontsize',12,'fontweight','bold');


%% Write new model
c11=vp.^2.*rho;
c55=vs.^2.*rho;

fid=fopen('c11_file','w');
fwrite(fid,c11','float32');
fclose(fid);

fid=fopen('c55_file','w');
fwrite(fid,c55','float32');
fclose(fid);


fid=fopen('rho_file','w');
fwrite(fid,rho','float32');
fclose(fid);


% fid=fopen('q_file','w');
% fwrite(fid,q','float32');
% fclose(fid);



%%
clearvars; close all; clc;

[Data,SuTraceHeaders,SuHeader]=ReadSu('vsp1.su','endian','b');
[m,n]=size(Data);
fclose('all');

dx=0.25;
dt=SuTraceHeaders(1).dt*1e-6;

t=(0:m-1)*dt;
%%
nn=100/dx;
z=(0:nn)*dx;
Ih=1    :nn;
Iv=Ih+n/2;
T=t<0.03;

figure;

subplot(121)
imagesc(t,z,tr_norm(Data(T,Ih),1)'); colormap(gray);
subplot(122)
imagesc(t,z,tr_norm(Data(T,Iv),1)'); colormap(gray);
