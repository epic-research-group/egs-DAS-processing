function [WAG,hel] = getWAG(R,alpha,L,ds,phi)

% function to construct G matrix (Chen Ning and Sava, 2018)
%  Usage: [WAG,hel] = getWAG(R,alpha,L,ds,phi)
%
%    Inputs: 
%         R- radius of mandrel (m)
%         alpha - wrapping angle of fiber (deg)
%         L - gauge length (m)
%         ds - distance between strain samples along fiber (only matters
%         for ouput array "hel")
%         phi - phase shift for helix (deg)
%
%    Outputs:
%         WAG  - WAG matrix for projecting strain onto fiber and averaging
%                   over gauge length
%         hel  - helix coordinates
%         

% first make helix with fiber length=L
dist=0:ds:L; % distance along fiber
t=dist/L;    % parameterization variable
h=L*sind(alpha); % length of axis
D=2*pi*R*tand(alpha); % axial distance of one wrap
ang=2*pi*h/D;         % angle for sines and cosines

phi=phi*pi/180;

hel=[R*cos(ang*t(:)+phi) R*sin(ang*t(:)+phi) h*t(:)];  % helix coordinates
%make sure straight fiber has x=0,y=0;
if alpha==90
    hel(:,1)=0;
end
%tangent vector's are
% tx=-R*ang*sin(ang*t(:)+phi);
% ty=R*ang*cos(ang*t(:)+phi);
% tz=h;

% magnitude of vectors:
EUC=sqrt(R^2*ang^2+h^2);

%G=[tx(:).^2 ty(:).^2 tz(:).^2 2*tx(:).*ty(:) 2*tx(:).*tz(:) 2*ty(:).*tz(:)]; 
% and WA just integrates from t= 0 to 1 along each column
% WAG
WAG = [R^2*ang*( ang*0.5 - (sin(2*(ang+phi)) - sin(2*phi) ) / (4));
       R^2*ang*( ang*0.5 + (sin(2*(ang+phi)) - sin(2*phi) ) / (4));
       h^2;
      -R^2*ang* ( sin(ang+phi)^2 - sin(phi)^2 );
       2*h*R* ( cos(ang+phi) - cos(phi) );
       2*h*R* ( sin(ang+phi) - sin(phi) )]' / EUC;

WAG=WAG/L; 
   





