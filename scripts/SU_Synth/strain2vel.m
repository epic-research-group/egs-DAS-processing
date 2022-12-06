function out=strain2vel(data,x,t,alpha,maxv)

% function to convert fiber strain-rate to particle velocity 
%   integrates strain-rate to get strain, then scales by phase velocity in
%   the f-k domain
%
% Usage:   out=fkspect(data,x,t)
%
% inputs:
%   data     = m x n data matrix one trace per column
%   x        = offset vector length(x) must equal n 
%   t        = time vector length(t) must equal m
%   
% returns:
%   out     = particle velocity


[m,n1]=size(data);

%padd in x
data=[zeros(m,2*n1) data zeros(m,2*n1)];
n=5*n1;

NFFT1 = (2^nextpow2(2^nextpow2(m))); 
NFFT2 = (2^nextpow2(2^nextpow2(n))); 

dt=t(2)-t(1);
dx=abs(mean(diff(x)));
df=1/(NFFT1*dt);
dk=1/(NFFT2*dx);
freq1 = (0:NFFT1-1)*df - (1/(2*dt));
kx = (0:NFFT2-1)*dk - (1/(2*dx));

% get phase velocities
[kx,freq]=meshgrid(kx,freq1); 
c= (freq.*kx) ./ (kx.^2 + 1e-8);
% c=freq./(kx);
% c(kx==0)=0;

if nargin>3
% figure out incident angles
inc=real(acosd(maxv./c));

% get sensitivity
fi=fiber_sens(alpha,inc,'p');
fi=abs(fi);



c=ifftshift(fi.*c);

else
    c=ifftshift(c);
end

%fft in time
%integrate first
data=cumtrapz(data)*dt;
fdata=fft2(data,NFFT1,NFFT2);

% now integrate in space
fdata=c.*fdata;

% ifft 
out=real(ifft2(fdata));
out=out(1:m,1:n);
st=2*n1+1;
en=st+n1-1;
out=out(:,st:en);








