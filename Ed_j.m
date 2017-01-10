clear all

% Variables :
k = physconst('Boltzman'); 
T = 300;
W = 4.5*10^(-6); %m
Wn_plus = 0.01e-6; %m
Wn = 0.9e-6; %m
e = 1.6*10^(-19);
Nd = 2*10^22; %m-3
Nd_plus = 8*10^24; %m-3
Lp_plus = 0.32e-6; %m    diffusion length of minority carriers in n+zone variation of minority-carrier diffusion length with carier concentration in GaAs liquid-phase epitaxial layers
Lp = 2e-6; %m            diffusion length of minority carriers in n zone
Ln =  20e-6; %m Dl of electron in n zone (problem no data  :s) 
epsilon = 12.9 *8.8541878176e-12 ; 
Dp_plus = 10e-4; %m2*s-1  depends on voltage ? which voltage
Dp = 10e-4; %m2*s-1 depends on voltage if any  ?
Sp = 3e3; % m*s-1
syms Vh x;
S = solve (e*Vh /(k*T) + log((Nd/Nd_plus) + 0.5*(e*Vh/(k*T))^2) == 0, Vh); %calculates potential across HLJ
Wa = sqrt(2*epsilon*k*T/(e^2*Nd)); % calculates vpa with previous Vh calculated
Senn_plus = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus));
Sen_plus_n = (Nd_plus/Nd)*(Dp/Lp)*coth(Wn/Ln);



b
