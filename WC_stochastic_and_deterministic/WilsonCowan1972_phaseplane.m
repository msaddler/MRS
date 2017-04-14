%WilsonCowan1972_.m
% variant with E is on the X-axis and I on the Y-axis
clear;
close all;
% constants
ke=1;ki=1;              % max value of response function p.10
%c1=12;c2=4;c3=13;c4=11; % coupling strengths between populations p. 10 Fig.4
c1=25;c2=4;c3=20;c4=2; % coupling strengths between populations Modified
%ae=1.2;Thetae=2.8; ai=1;Thetai=4;    % Boltzman function parms Fig.4
ae=.8;Thetae=2.9; ai=5;Thetai=1;    % Boltzman function parms Modified
re=1;ri=1;              % refractory period
%P=0; Q=0;               % external input
P=0.1; Q=0;               % external input Modified
I=-0.1:0.01:0.5;            % ranges for I (inhibitory firing rate)and
E=I;                        % excitatory firing rate E
[Ig,Eg]=meshgrid(I,E);  % this creates the grid for the quiver command
lI=length(I);
lE=length(E);
for i=1:lI;
    for j=1:lE;
        help=c1*E(j) - c2*I(i) + P;
        help=1/(1 + exp(-ae*(help-Thetae)));
        help=help-(1/(1+exp(ae*Thetae)));   %Eq (15), p.11
        dE(i,j)=-E(j) + (ke-re*E(j)) * help;
       
        help=c3*E(j) - c4*I(i) + Q;
        help=1/(1 + exp(-ai*(help-Thetai)));
        help=help-(1/(1+exp(ai*Thetai)));   %Eq (15), p.11
        dI(i,j)=-I(i) + (ki-ri*I(i)) * help;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%
% plot results
% null-clines
% We find the n-null cline
% using the contour command
figure;
%contour(I,E,dE',[0 0],'g');
contour(E,I,dE,[0 0],'g');
hold;
%contour(I,E,dI',[0 0],'r')
contour(E,I,dI,[0 0],'r')
% make the mesh sparser for display purposes
step=2;
Es=Eg(1:step:lE,1:step:lI);
Is=Ig(1:step:lE,1:step:lI);
dEs=dE(1:step:lE,1:step:lI);
dIs=dI(1:step:lE,1:step:lI);
%quiver(Is,Es,dIs',dEs',3); 
quiver(Es,Is,dEs',dIs',3); 
 
title('Fig. 4 Wilson and Cowan 1972 green - E-null; red - I-null ')
%xlabel('I')
%ylabel('E')
xlabel('E')
ylabel('I')