clear all

% Variables :
k = physconst('Boltzman'); 
T = 333;
W = 4.5*10^(-6);
e = 1.6*10^(-19);
Nd = 8*10^22; %m-3
Nd_plus = 2*10^23; %m-3
epsilon = 12.9 ; 
syms Vh x;
S = solve (e*Vh /(k*T) + log((Nd/Nd_plus) + 0.5*(e*Vh/(k*T))^2) == 0, Vh); %calculates potential across HLJ
Wa = sqrt(2*epsilon*k*T/(e^2*Nd))*atan(e*S/(sqrt(2)*k*T)*sqrt(Nd_plus/Nd)); % calculates vpa with previous Vh calculated







