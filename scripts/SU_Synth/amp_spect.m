function [spect,freq,an]=amp_spect(data,dt)


[m,~]=size(data);

%fft stuff
nfft=2^nextpow2(m);
fs=1/dt;
nyq=1/(2*dt);
df=fs/nfft;
freq1=((0:nfft-1))*df-nyq;
DD=fft(data,nfft,1);
DD=fftshift(DD,1);
spect=2*(abs(DD)/m);
an=angle(DD);

I=freq1>=0;
spect=spect(I,:);
freq=freq1(I);

