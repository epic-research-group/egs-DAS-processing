%%%
%%%  Iterative least squares algorthim to find L1 norm 

% This is the IRLS algorithm described on p 34 of Dick Asters Inverse theory book
% A=matrix (physics of your problem)
% d=data vector
% T=tolerance (Tau in book) 
% e=minimum residual (eta in book)
% flag will force a zero intercept if it is equal to the column of A that
% represents the intercept value.  Default is 0.

function [c,d,residual] = L1(A,d,T,e,maxiter)

if nargin<4
    error('L1 wants at least four arguments');
end

if nargin == 4;
   
    maxiter = 100;
    fl = 0;

end

if nargin == 5;
    fl = 0;
end
fl = 0;


delta=10^10;

    if fl;
        m2 = lsqr(A,d,1e-2);
    else
        m2=A\d;
    end
    r=d-A*m2;
    r(r==0) = e; %R(1:length(D)) = D;
    R=diag(abs(r.^-1));
    m1=m2;

cnt = 0;
while delta>T;
    
    cnt = cnt +1;
      
    if cnt > maxiter;
        
        c=m1;
        d=A*m1;
        residual=r;
        return
        
    end
    
    m_old=m1;
    ATR=A'*R;
    
    if fl;
        m1 = lsqr(ATR*A,ATR*d,1e-2);
    else
        m1=(ATR*A)\(ATR*d); 
    end
        r=d-A*m1; %R(1:length(D)) = D;
        I=find(r==0);
        if isempty(I)==0;
            r(I)=e;
        end
            
    
    R=diag(abs(r.^-1));
    
    a=(m1-m_old)'*(m1-m_old);
    b=m1'*m1;
    
    delta=a/(1+b);
    
    
    
end

c=m1;
d=A*m1;
residual=r;

