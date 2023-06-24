function C = tr_norm(A,dir)

if exist('dir','var') == 0
    
    dir = 1;
    
end


% normalize seismic traces

C = zeros(size(A));

if dir == 1

    for i = 1:length(A(1,:))
    

        C(:,i) = A(:,i)/max(abs(A(:,i)));

    
    end
    
elseif dir == 2
    
    A = A';
    C = zeros(size(A));
    
    for i = 1:length(A(1,:))
    
    
        C(:,i) = A(:,i)/max(abs(A(:,i)));

    end
    
    C = C';
    
    elseif dir==3
        C=A/max(abs(A(:)));
end


C(isnan(C))=0;


    

