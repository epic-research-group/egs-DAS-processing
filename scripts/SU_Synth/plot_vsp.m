clear all; close all; clc

% load su data
[Data1,SuTraceHeaders,SuHeader]=ReadSu('./vsp.su');%,'endian','b');
[Data,SuTraceHeaders,SuHeader]=ReadSu('./vsp1.su');%,'endian','b');
[m,n]=size(Data);
fclose('all'); % for some reason ReadSu doesn't close the files

dx=0.25; % x spacing between vsp.su and vsp1.su
dz=0.25; % z spacing between synthetic traces
dt=SuTraceHeaders(1).dt*1e-6;

t=(0:m-1)*dt;

col=1200;
row=1200;

% get the upper 100 meters
x=(0:400)*dz;
Ih=1:length(x);
Iv=Ih+col;
Datah1=Data1(:,Ih); Datah=Data(:,Ih);
Datav1=Data1(:,Iv); Datav=Data(:,Iv);
Datavs=Datav1; % save original velocity data
Datahs=Datah1;
%%
% get amplitude spectra
[sph,f]=amp_spect(Datah,dt);
[spv,f]=amp_spect(Datav,dt);

% plot data and amplitude spectra
figure;
subplot(221);
imagesc(t,x,(Datah)');
xlim([0 0.02])

colormap(gray);
subplot(222);
imagesc(f,x,log10(sph)');
% xlim([0 2e4]);

subplot(223);
imagesc(t,x,(Datav)');
xlim([0 0.02])

subplot(224);
imagesc(f,x,log10(spv)');
% xlim([0 1e4]);

%% Let's get strain
% first let's go from velocity to displacement
Datav1=cumsum(Datav1)*dt;
Datav=cumsum(Datav)*dt;
Datah1=cumsum(Datah1)*dt;
Datah=cumsum(Datah)*dt;

svz=[0*Datav1(:,1) diff(Datav1,1,2)/dz]; % du_vertical/dz
shx=[Datah1-Datah]/dx; % du_horizontal/dx

svx=[Datav1-Datav]/dx; % du_vertical/dx
shz=[0*Datah1(:,1) diff(Datah1,1,2)/dz]; % du_horizontal/dz

%e11=shx;
%e33=svz;
%e13=e31=0.5*(svx+shz)

%% Project strain onto 35 degree helical fiber
R=0.0125; % radius of mandrel (m)
alpha=35; % wrapping angle
L=2;      % gauge length
ds=0.001; % sampling interval (only matters for plotting helix)
phi=0;    % phase shift for helix

% get WAG to project strain tensor onto fiber
[WAG,hel] = getWAG(R,alpha,L,ds,phi);

%make a new data set
xnew=1:1:100;
SS=zeros(length(t),length(xnew)); % strain

for i=1:length(t)
    
    for j=1:length(xnew)
        
        I=abs(x-xnew(j))<=L/2; % channels within one gauge length of center
        n=sum(I);
        %strain tensor
        s=zeros(6,n);
        s(1,:)=shx(i,I);
        s(3,:)=svz(i,I);
        s(5,:)=0.5*(svx(i,I) +shz(i,I));
        SS(i,j)=mean(WAG*s,2); % average response
        
        
    end
    
end

% plot the 35 degree to show no shear waves
figure;
subplot(131);
imagesc(xnew,t,tr_norm(Datavs)); colormap(gray);
xlabel('meters');
ylabel('seconds');
title('Vertical Component velocity');

subplot(132);
imagesc(xnew,t,tr_norm(Datahs)); colormap(gray);
xlabel('meters');
ylabel('seconds');
title('Horizontal Component velocity');

subplot(133);
imagesc(xnew,t,tr_norm(SS)); colormap(gray);
xlabel('meters');
ylabel('seconds');
title('35 deg helical fiber');

%%
R=0.0125; % radius of mandrel (m)
alpha=60; % wrapping angle (90 = straight)
L=2;      % gauge length
ds=0.001; % sampling interval (only matters for plotting helix)
phi=0;    % phase shift for helix

% get WAG to project strain tensor onto fiber
[WAG,hel] = getWAG(R,alpha,L,ds,phi);

%now make a new data set
xnew=1:1:100;

SS=zeros(length(t),length(xnew)); % strain
HH=SS; %Horizontal velocity
VV=SS; %Vertical velocity
for i=1:length(t)
    
    
    for j=1:length(xnew)
        
        I=abs(x-xnew(j))<=L/2;
        n=sum(I);
        %strain tensor
        s=zeros(6,n);
        mind=round(n/2); % index of geophone in center of gauge
        s(1,:)=shx(i,I);
        s(3,:)=svz(i,I);
        s(5,:)=0.5*(svx(i,I) +shz(i,I));
        tmp=Datahs(i,I);
        HH(i,j)=tmp(mind);
        tmp=Datavs(i,I);
        VV(i,j)=tmp(mind);
        SS(i,j)=mean(WAG*s,2);
    
    end
    
end


%% try to recover velocity

out1=strain2vel(diff(SS)/dt,xnew,t,90,8000); % f-k version
out=strain2vel3(diff(SS)/dt,2,1,0.2);       % inverse of difference operator

%%
figure;
h1=subplot(141);
imagesc(xnew,t,tr_norm(SS,3)); colormap(gray);
title('Strain Rate');

h2=subplot(142);
imagesc(xnew,t,tr_norm(VV));
title('Vertical particle velocity');

h3=subplot(143);
imagesc(xnew,t,tr_norm(out));
title('Inverse problem');

h4=subplot(144);
imagesc(xnew,t,tr_norm(out1));
title('F-k');
% colorbar
linkaxes([h1 h2 h3 h4]);
%%
% compare results to true velocity
figure;
for i=1:length(xnew)
%     subplot(211);
    plot(t,(VV(:,i)),'b')
    hold on;
    plot(t(2:end),(out(:,i)),'r');
    plot(t(2:end),(out1(:,i)),'g');
    xlim([0 0.025])
    pause;
    hold off;
end

%%








