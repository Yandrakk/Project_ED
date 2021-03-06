clear all

% Variables :
k = physconst('Boltzman'); 
T = 300;
W = 4.5*10^(-6); %m
Wn_plus = 0.01e-6; %m
Wn = 0.9e-6; %m
e = 1.6*10^(-19);
range_ndplus = [10^22 : 10^20 : 10^26];
Nd = 2e22*ones(1,length(range_ndplus)); %m-3
% Nd_plus = 8*10^24; %m-3

Lp_plus = 0.32e-6; %m    diffusion length of minority carriers in n+zone variation of minority-carrier diffusion length with carier concentration in GaAs liquid-phase epitaxial layers
Lp = 2e-6; %m            diffusion length of minority carriers in n zone
Ln =  20e-6; %m Dl of electron in n zone (problem no data  :s) 
epsilon = 12.9 *8.8541878176e-12 ; 
Dp_plus = 400e-4 * k*T/e; %� * k*T./e
Dp = Dp_plus; % most likely we don't care about doping 
Sp = 3e3; % m*s-1
syms Vh x, 
<<<<<<< HEAD
%S = solve (e*Vh /(k*T) + log((Nd./range_ndplus) + 0.5*(e*Vh/(k*T))^2) == 0, Vh); %calculates potential across HLJ
%Wa = sqrt(2*epsilon*k*T./(e^2*Nd))*atan(e*S/(sqrt(2)*k*T)*sqrt(Nd_plus/Nd)); % calculates vpa with previous Vh 

S = solve (e*Vh /(k*T) + log((Nd./range_ndplus) + 0.5*(e*Vh/(k*T))^2) == 0, Vh); %calculates potential across HLJ
Wa = sqrt(2*epsilon*k*T./(e^2*Nd))*atan(e*S/(sqrt(2)*k*T)*sqrt(range_ndplus./Nd)); % calculates vpa with previous Vh
%calculated NEED TO ADD SIN
>>>>>>> origin/master
Senn_plus = (Nd./range_ndplus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus));
Sen_plus_n = (range_ndplus./Nd)*(Dp/Lp)*coth(Wn/Ln);
Fh_1 = 1./(1+(Senn_plus./Sen_plus_n).*(range_ndplus./Nd));
Fh_1_corrected= Fh_1 *sech(Wn/Lp);
figure;
semilogx(range_ndplus,Fh_1,range_ndplus,Fh_1_corrected);

 





