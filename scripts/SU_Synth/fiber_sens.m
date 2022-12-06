function out=fiber_sens(alpha,th,type)

%function to predict fiber sensitivity as a function of incidence angle and
%wrapping angle
% 
% Inputs:
%       alpha - wrapping angle (90=straight) (deg)
%       th    - incidence angle (deg)
%       type  = 'p' or 's'



th=th*pi/180-pi/2;

% rotation matrix
R11=cos(th);
R13=sin(th);
R22=1;
R31=-R13;
R33=R11;


%         %R*e
%         tmp11=R11*e11 + R13*e31;
%         tmp21=R22*e21;
%         tmp31=R31*e11 + R33*e31;
% 
%         tmp12=R11*e12 + R13*e32;
%         tmp22=R22*e22;
%         tmp32=R31*e12 + R33*e32;
% 
%         tmp13=R11*e13 + R13*e33;
%         tmp23=R22*e23;
%         tmp33=R31*e13 + R33*e33;
% 
%         % R*e*R'
%         e11=tmp11*R11 + tmp13*R13;
%         e22=tmp22*R22;
%         e33=tmp31*R31 + tmp33*R33;


switch type
    case 'p'
%         e1=[0 0 0; 0 0 0; 0 0 1]; % pwave
        %R*e
        e33p=1;
        tmp13=R13*e33p;
        tmp33=R33*e33p;

        % R*e*R'
        e11=tmp13.*R13;
        e22=0;
        e33=tmp33.*R33;
        
        
    case 's'
%         e1=[0 0 1; 0 0 0; 1 0 0]; % swave
        e13s=1;
        e31s=1;
        %R*e
        tmp11=R13*e31s;
        tmp31=R33*e31s;
        
        tmp13=R11*e13s;
        tmp33=R31*e13s;

        % R*e*R'
        e11=tmp11.*R11 + tmp13.*R13;
        e22=0;
        e33=tmp31.*R31 + tmp33.*R33;
end




    out=e11*sind(alpha)^2 + 0.5*e22*cosd(alpha)^2 + ...
        0.5*e33*cosd(alpha)^2;

end