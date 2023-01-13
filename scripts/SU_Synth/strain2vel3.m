function out = strain2vel3(in,Lg,dx,dmp)

% function to do an inverse problem to get velocity from strain-rate
% in - input (strain-rate)
% Lg - gauge length
% dx - channel interval
% dmp - damping paramter

[m,n]=size(in);
x=(0:n-1)*dx; % coordinates for strain-rate
newx=(-Lg/2:dx:(max(x)+Lg/2));
ind=Lg/dx;
L=zeros(n,length(newx));

for i=1:n
    
    L(i,i)=-1/Lg;
    L(i,i+ind)=1/Lg;
    
end

out=zeros(size(in));
LHS=[L; dmp*diag(ones(length(newx),1))];
for i=1:m
    
    RHS=[in(i,:)'; zeros(length(newx),1)];
    [m,~,~]=L1(LHS,RHS,5e-6,1e-2,100);
    
    out(i,:)=interp1(newx,m',x);
    
end