clear all

% Variables :
k = physconst('Boltzman'); 
T = 300;
W = 4.5*10^(-6); %m
Wn_plus = 0.01e-6; %m
Wn = 0.9e-6; %m
Wp = 3e-6; %m
Wn_range = [0.001e-6 : 5e-9 : 5e-6]; %wn range pour tracer photocurrent density en fonction de Wn
e = 1.6*10^(-19);
lambda = [0.24: 0.01 : 1] ; % en �
lambda_fixe = 0.4; % en � on prend celle qui donne le plus de courant 
range_ndplus = [10^22 : 10^20 : 10^26];
Nd_plus = 2e25 ; % m-3 
Nd = 2e22; %m-3
Na = 2e25; %m-3
tau_p = 0.1e-6/sqrt((Nd*10^(-6)-10^16)/10^16) ; %formule photo
tau_p_plus = 0.1e-6/sqrt((Nd_plus*10^(-6)-10^16)/10^16) ; %formule photo
Dp = 26e-7/(2.5e-3 + 4e-21*Nd*10^(-6)); %photo formula
Dp_plus = 26e-7/(2.5e-3 + 4e-21*Nd_plus*10^(-6));    
Lp_plus = sqrt(Dp_plus*tau_p_plus); %m    diffusion length of minority carriers in n+zone variation of minority-carrier diffusion length with carier concentration in GaAs liquid-phase epitaxial layers
Lp = sqrt(Dp*tau_p); %m            diffusion length of minority carriers in n zone
tau_n = 0.1e-6/sqrt((Na*10^(-6)-10^16)/10^16);
Dn = 200e-4; % m2.s-1 electron diffusivit� (in p z  one )
Ln =  sqrt(Dn*tau_n); 
epsilon = 12.9 *8.8541878176e-12 ; %relative permittivity gaas

  

Sp = 3e3; % m*s-1
Sp_plus = 4e3; 
Sb = 1e5 ;  % surface recombination velocity in p-type gaas
syms Vh ; 
S = solve (e*Vh /(k*T) + log((Nd/Nd_plus) + 0.5*(e*Vh/(k*T))^2) == 0, Vh); %calculates potential across HLJ
Wa = sqrt(2*epsilon*k*T./(e^2*Nd))*atan(e*S/(sqrt(2)*k*T)*sqrt(Nd_plus/Nd)); % calculates vpa with previous Vh

Senn_plus = (Nd./range_ndplus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus));
Sen_plus_n = (range_ndplus./Nd)*(Dp/Lp)*coth(Wn/Ln);
Fh_1 = 1./(1+(Senn_plus./Sen_plus_n).*(range_ndplus./Nd));
Fh_1_corrected= Fh_1 *sech(Wn/Lp);
%figure;
%semilogx(range_ndplus,Fh_1,range_ndplus,Fh_1_corrected);
f_lambda = arrayfun(@f, lambda); 
a_lambda = arrayfun(@a, lambda);
 
% Photocurrent density en fonction de lambda   
Senn_plus_lambda = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %ne depend pas de lambda donc blc
Sen_plus_n_lambda = (Nd_plus/Nd)*(Dp/Lp)*coth(Wn/Ln); %depend pas de lambda non plus
Fh_1_lambda = 1/(1+(Senn_plus_lambda./Sen_plus_n_lambda)*(Nd_plus/Nd));
Fh_1_corrected_lambda= Fh_1_lambda *sech(Wn/Lp); 

Tau = f_lambda.*((1-0.05)*e*Lp_plus*a_lambda)./(-1+(a_lambda*Lp_plus).^2); % R = 0.05, black material, as your soul
Jn0_plus_x1_lambda = Fh_1_corrected_lambda*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus))));
J_Wa_x3_lambda = e*f_lambda.*(1-0.05).*exp(-a_lambda*Wn_plus).*(1-exp(-a_lambda*Wa)).*sech(Wn/Lp);
k1_lambda = e * f_lambda .* (1-0.05) .* a_lambda.* Lp./((a_lambda*Lp).^2 -1 ); 
k2_lambda = (e*f_lambda.*(1-0.05).*a_lambda*Ln./(((a_lambda*Ln).^2) -1)).*exp(-a_lambda*Wp);
J_p_lambda = -k1_lambda.*Lp.*a_lambda.*exp(-a_lambda.*Wn) + (k1_lambda./((Senn_plus_lambda.*Lp/Dp).*sinh(Wn/Lp)+cosh(Wn/Lp))).*((Senn_plus_lambda.*Lp/Dp)+ a_lambda.*Lp - exp(-a_lambda.*Wn).*((Senn_plus_lambda*Lp/Dp).*cosh(Wn/Lp) + sinh(Wn/Lp))); % photocurrent generated by p layer calculated at x3
J_n_lambda = k2_lambda .* a_lambda .* Ln - (k2_lambda./((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln))).*((Sb*Ln/Dn).*(cosh(Wp/Ln) - exp(-a_lambda.*Wp)) + sinh(Wp/Ln) + a_lambda.*Ln.*exp(-a_lambda*Wp)) ;

J_Wd_lambda = e .* f_lambda.*(1-0.05).*exp(-a_lambda*Wn).*(1-exp(-a_lambda*W));
J_tot_lambda = Jn0_plus_x1_lambda + J_Wa_x3_lambda +J_p_lambda +J_n_lambda +J_Wd_lambda;
figure;
plot(lambda, J_tot_lambda, 'md');
figure;
plot(lambda, J_p_lambda,'bo', lambda, J_n_lambda, 'gx',lambda, J_Wd_lambda, 'rs',lambda, Jn0_plus_x1_lambda,'yv', lambda,J_Wa_x3_lambda,'kd' );

%Photocurent density en fonction de Wn
Senn_plus_Wn = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %ne depend pas de Wn donc blc
Sen_plus_n_Wn = (Nd_plus/Nd)*(Dp/Lp)*coth(Wn_range/Ln); %depend de Wn donc matrice
Fh_1_Wn = 1./(1+(Senn_plus_Wn./Sen_plus_n_Wn).*(Nd_plus/Nd));
Fh_1_corrected_Wn= Fh_1_Wn .*sech(Wn_range./Lp); %matrice
Tau_lambda_fixe = f(lambda_fixe)*((1-0.05)*e*Lp_plus*a(lambda_fixe))/(-1+(a(lambda_fixe)*Lp_plus)^2); 
Jn0_plus_x1_Wn = Fh_1_corrected_Wn.*(Tau_lambda_fixe * ((-a(lambda_fixe)*Lp_plus).*exp(-a(lambda_fixe)*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a(lambda_fixe)*Lp_plus - exp(-a(lambda_fixe)*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))));
J_Wa_x3_Wn = e*f(lambda_fixe)*(1-0.05)*exp(-a(lambda_fixe)*Wn_plus)*(1-exp(-a(lambda_fixe)*Wa))*sech(Wn_range/Lp);

k1 = e * f(lambda_fixe) * (1-0.05) * a(lambda_fixe)* Lp/((a(lambda_fixe)*Lp)^2 -1 ); 
k2 = (e*f(lambda_fixe)*(1-0.05)*a(lambda_fixe)*Ln/(((a(lambda_fixe)*Ln)^2) -1))*exp(-a(lambda_fixe)*Wp);
J_p_Wn = -k1.*Lp.*a(lambda_fixe).*exp(-a(lambda_fixe).*Wn_range) + (k1./((Senn_plus_Wn.*Lp/Dp).*sinh(Wn_range/Lp)+cosh(Wn_range/Lp))).*((Senn_plus_Wn.*Lp/Dp)+ a(lambda_fixe).*Lp - exp(-a(lambda_fixe).*Wn_range).*((Senn_plus_Wn*Lp/Dp).*cosh(Wn_range/Lp) + sinh(Wn_range/Lp))); % photocurrent generated by p layer calculated at x3
J_n_Wn = k2 * a(lambda_fixe) * Ln - (k2/((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln)))*((Sb*Ln/Dn)*(cosh(Wp/Ln) - exp(-a(lambda_fixe)*Wp)) + sinh(Wp/Ln) + a(lambda_fixe)*Ln*exp(-a(lambda_fixe)*Wp)) ;

J_Wd = e * f(lambda_fixe)*(1-0.05)*exp(-a(lambda_fixe)*Wn_range)*(1-exp(-a(lambda_fixe)*W));
J_tot = J_p_Wn + J_n_Wn  + Jn0_plus_x1_Wn + J_Wa_x3_Wn + J_Wd ;
figure;
plot(Wn_range, J_p_Wn,'bo', Wn_range, J_n_Wn, 'gx',Wn_range, J_Wd, 'rs',Wn_range, Jn0_plus_x1_Wn,'yv', Wn_range,J_Wa_x3_Wn,'kd' );
xlabel('Wn (m)');
ylabel ('Photocurrent mA/cm�');
figure;
plot(Wn_range, J_tot,'mp');
xlabel('Wn (m)');
ylabel ('Photocurrent mA/cm�');

%Internal quantum efficiency VS lambda

 
 











