clear all
close all

% Variables :
k = physconst('Boltzman'); 
h = 6.63e-34; 
c = physconst('lightspeed');
T = 300;
epsilon = 12.9 *8.8541878176e-12 ; %relative permittivity gaas

W = 4.5*10^(-6); %m
Wn_plus = 0.01e-6; %m
Wn_plus_1 = 0.02e-6;
Wn_plus_2 = 0.03e-6;
Wn_plus_3 = 0.04e-6;
Wn_plus_4 = 0.05e-6;
Wn_plus_range = [0.001e-6 : 0.01e-9 : 0.2e-6];
Wn = 0.9e-6; %m
Wp = 2e-6; %m
Wn_range = [0.001e-6 : 5e-9 : 5e-6]; %wn range pour tracer photocurrent density en fonction de Wn
e = 1.6*10^(-19);
q=e;
lambda = [0.26: 0.01 : 1] ; % en µ
lambda_fixe = 0.7; % en µ on prend celle qui donne le plus de courant 
range_ndplus = [10^22 : 10^20 : 10^26];
Nd_plus = 5e25 ;
Nd_plus_1 = 2e23;
Nd_plus_2 = 2e22;
% m-3 
Nd = 2e22; %m-3
Nd_1 = 6e23;
Nd_2 = 1e25;
Na = 2e25; %m-3
tau_p = 0.1e-6/sqrt((Nd*10^(-6)-10^16)/10^16) ; %formule photo
tau_p_plus = 0.1e-6/sqrt((Nd_plus*10^(-6)-10^16)/10^16) ; %formule photo
Dp = 26e-7/(2.5e-3 + 4e-21*Nd*10^(-6)); %photo formula
Dp_plus = 26e-7/(2.5e-3 + 4e-21*Nd_plus*10^(-6));    
Lp_plus = sqrt(Dp_plus*tau_p_plus); %m    diffusion length of minority carriers in n+zone variation of minority-carrier diffusion length with carier concentration in GaAs liquid-phase epitaxial layers
Lp = sqrt(Dp*tau_p); %m            diffusion length of minority carriers in n zone
tau_n = 0.1e-6/sqrt((Na*10^(-6)-10^16)/10^16);
Dn = 200e-4; % m2.s-1 electron diffusivité (in p z  one )
Ln =  sqrt(Dn*tau_n);
ni = 1.8e12; %m-3
V_bi = (k*T/q)*log(Na/ni) + (k*T/q)*log(Nd/ni);
Wd = sqrt((2*epsilon/q)*((Na+Nd)/(Na*Nd))*(V_bi));
  

Sp = 3e3; % m*s-1
Sp_plus = 4e3; 
Sb = 1e6 ;  % surface recombination velocity in p-type gaas
syms Vh;

S = solve (e*Vh /(k*T) + log((Nd/Nd_plus) + 0.5*(e*Vh/(k*T))^2) == 0, Vh); %calculates potential across HLJ
Wa = sqrt(2*epsilon*k*T./(e^2*Nd))*atan(e*S/(sqrt(2)*k*T)*sqrt(Nd_plus/Nd)); % calculates vpa with previous Vh

Senn_plus = (Nd./range_ndplus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus));
Sen_plus_n_Nd_plus = (range_ndplus./Nd)*(Dp/Lp)*coth(Wn/Ln);
Fh_1 = 1./(1+(Senn_plus./Sen_plus_n_Nd_plus).*(range_ndplus./Nd));
Fh_1_corrected= Fh_1 *sech(Wn/Lp);
Senn_plus_Wn_plus = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus_range / Lp_plus))./(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus_range/Lp_plus)); %effect of n+ on the n carriers 
Sen_plus_n_Wn_plus = (Nd_plus./Nd)*(Dp/Lp)*coth(Wn/Ln);
%figure;
%plot(Wn_plus_range, Senn_plus_Wn_plus);
Fh_1_Wn_plus = 1./(1+(Senn_plus_Wn_plus./Sen_plus_n_Wn_plus).*(Nd_plus./Nd));
Fh_1_Wn_plus_corrected= Fh_1_Wn_plus *sech(Wn/Lp);
figure;
plot(Wn_plus_range, Fh_1_Wn_plus,'r-', Wn_plus_range, Fh_1_Wn_plus_corrected,'g-');
title('High Low Junction Factor VS Wn+');
xlabel('Wn+ (m)');
ylabel('Fh_1');
legend('Fh_1', 'Fh_1 corrected');
figure;
semilogx(range_ndplus,Fh_1,range_ndplus,Fh_1_corrected);
title('High Low Junction Factor VS Nd+ doping');
xlabel('Nd+ (m-3)');
ylabel('Fh_1');


legend('Fh_1', 'Fh_1 corrected');
f_lambda = arrayfun(@f, lambda); 
a_lambda = arrayfun(@a, lambda);
R_lambda = arrayfun(@R, lambda);

 
% Photocurrent density vs lambda
Senn_plus_lambda = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %effect of n+ on the n carriers 
Senn_plus_lambda_1 = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus_1 / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus_1/Lp_plus));
Senn_plus_lambda_2 = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus_2 / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus_2/Lp_plus));
Senn_plus_lambda_3 = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus_3 / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus_3/Lp_plus));
Senn_plus_lambda_4 = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus_4 / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus_4/Lp_plus));
Sen_plus_n_lambda = (Nd_plus/Nd)*(Dp/Lp)*coth(Wn/Ln); %effect of the n region on the n+ carriers 
Fh_1_lambda = 1/(1+(Senn_plus_lambda./Sen_plus_n_lambda)*(Nd_plus/Nd));
Fh_1_lambda_1 = 1/(1+(Senn_plus_lambda_1./Sen_plus_n_lambda)*(Nd_plus/Nd));
Fh_1_lambda_2 = 1/(1+(Senn_plus_lambda_2./Sen_plus_n_lambda)*(Nd_plus/Nd));
Fh_1_lambda_3 = 1/(1+(Senn_plus_lambda_3./Sen_plus_n_lambda)*(Nd_plus/Nd));
Fh_1_lambda_4 = 1/(1+(Senn_plus_lambda_4./Sen_plus_n_lambda)*(Nd_plus/Nd));
Fh_1_corrected_lambda= Fh_1_lambda *sech(Wn/Lp); 
Fh_1_corrected_lambda_1= Fh_1_lambda_1 *sech(Wn/Lp); 
Fh_1_corrected_lambda_2= Fh_1_lambda_2 *sech(Wn/Lp); 
Fh_1_corrected_lambda_3= Fh_1_lambda_3 *sech(Wn/Lp); 
Fh_1_corrected_lambda_4= Fh_1_lambda_4 *sech(Wn/Lp); 

Tau = f_lambda.*((1-R_lambda).*e.*Lp_plus.*a_lambda)./(-1+(a_lambda.*Lp_plus).^2); % R = R, black material
Jn0_plus_x1_lambda = Fh_1_corrected_lambda*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))); 
Jn0_plus_x1_lambda_1 = Fh_1_corrected_lambda_1*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus_1) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus_1)*((Sp*Lp_plus*cosh(Wn_plus_1/Lp_plus)/Dp_plus) + sinh(Wn_plus_1/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus_1/Lp_plus)/Dp_plus) + cosh(Wn_plus_1/Lp_plus))));
Jn0_plus_x1_lambda_2 = Fh_1_corrected_lambda_2*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus_2) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus_2)*((Sp*Lp_plus*cosh(Wn_plus_2/Lp_plus)/Dp_plus) + sinh(Wn_plus_2/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus_2/Lp_plus)/Dp_plus) + cosh(Wn_plus_2/Lp_plus))));
Jn0_plus_x1_lambda_3 = Fh_1_corrected_lambda_3*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus_3) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus_3)*((Sp*Lp_plus*cosh(Wn_plus_3/Lp_plus)/Dp_plus) + sinh(Wn_plus_3/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus_3/Lp_plus)/Dp_plus) + cosh(Wn_plus_3/Lp_plus))));
Jn0_plus_x1_lambda_4 = Fh_1_corrected_lambda_4*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus_4) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus_4)*((Sp*Lp_plus*cosh(Wn_plus_4/Lp_plus)/Dp_plus) + sinh(Wn_plus_4/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus_4/Lp_plus)/Dp_plus) + cosh(Wn_plus_4/Lp_plus))));
J_Wa_x3_lambda = e*f_lambda.*(1-R_lambda).*exp(-a_lambda*Wn_plus).*(1-exp(-a_lambda*vpa(Wa))).*sech(Wn/Lp);
J_Wa_x3_lambda_1 = e*f_lambda.*(1-R_lambda).*exp(-a_lambda*Wn_plus_1).*(1-exp(-a_lambda*vpa(Wa))).*sech(Wn/Lp);
J_Wa_x3_lambda_2 = e*f_lambda.*(1-R_lambda).*exp(-a_lambda*Wn_plus_2).*(1-exp(-a_lambda*vpa(Wa))).*sech(Wn/Lp);
J_Wa_x3_lambda_3 = e*f_lambda.*(1-R_lambda).*exp(-a_lambda*Wn_plus_3).*(1-exp(-a_lambda*vpa(Wa))).*sech(Wn/Lp);
J_Wa_x3_lambda_4 = e*f_lambda.*(1-R_lambda).*exp(-a_lambda*Wn_plus_4).*(1-exp(-a_lambda*vpa(Wa))).*sech(Wn/Lp);%photocurrent of space charge zone (n-n+ zone)
k1_lambda = e * f_lambda .* (1-R_lambda) .* a_lambda.* Lp./((a_lambda*Lp).^2 -1 ); 
k2_lambda = (e*f_lambda.*(1-R_lambda).*a_lambda*Ln./(((a_lambda*Ln).^2) - 1)).*exp(-a_lambda*(Wn_plus+Wn+Wa+Wd));
k2_lambda_1 = (e*f_lambda.*(1-R_lambda).*a_lambda*Ln./(((a_lambda*Ln).^2) - 1)).*exp(-a_lambda*(Wn_plus_1+Wn+Wa+Wd));
k2_lambda_2 = (e*f_lambda.*(1-R_lambda).*a_lambda*Ln./(((a_lambda*Ln).^2) - 1)).*exp(-a_lambda*(Wn_plus_2+Wn+Wa+Wd));
k2_lambda_3 = (e*f_lambda.*(1-R_lambda).*a_lambda*Ln./(((a_lambda*Ln).^2) - 1)).*exp(-a_lambda*(Wn_plus_3+Wn+Wa+Wd));
k2_lambda_4 = (e*f_lambda.*(1-R_lambda).*a_lambda*Ln./(((a_lambda*Ln).^2) - 1)).*exp(-a_lambda*(Wn_plus_4+Wn+Wa+Wd));

J_n_lambda = -k1_lambda.*Lp.*a_lambda.*exp(-a_lambda.*Wn) + (k1_lambda./((Senn_plus_lambda.*Lp/Dp).*sinh(Wn/Lp)+cosh(Wn/Lp))).*((Senn_plus_lambda.*Lp/Dp)+ a_lambda.*Lp - exp(-a_lambda.*Wn).*((Senn_plus_lambda*Lp/Dp).*cosh(Wn/Lp) + sinh(Wn/Lp))); % photocurrent generated by n layer
J_n_lambda_1 = -k1_lambda.*Lp.*a_lambda.*exp(-a_lambda.*Wn) + (k1_lambda./((Senn_plus_lambda_1.*Lp/Dp).*sinh(Wn/Lp)+cosh(Wn/Lp))).*((Senn_plus_lambda_1.*Lp/Dp)+ a_lambda.*Lp - exp(-a_lambda.*Wn).*((Senn_plus_lambda_1*Lp/Dp).*cosh(Wn/Lp) + sinh(Wn/Lp)));
J_n_lambda_2 = -k1_lambda.*Lp.*a_lambda.*exp(-a_lambda.*Wn) + (k1_lambda./((Senn_plus_lambda_2.*Lp/Dp).*sinh(Wn/Lp)+cosh(Wn/Lp))).*((Senn_plus_lambda_2.*Lp/Dp)+ a_lambda.*Lp - exp(-a_lambda.*Wn).*((Senn_plus_lambda_2*Lp/Dp).*cosh(Wn/Lp) + sinh(Wn/Lp)));
J_n_lambda_3 = -k1_lambda.*Lp.*a_lambda.*exp(-a_lambda.*Wn) + (k1_lambda./((Senn_plus_lambda_3.*Lp/Dp).*sinh(Wn/Lp)+cosh(Wn/Lp))).*((Senn_plus_lambda_3.*Lp/Dp)+ a_lambda.*Lp - exp(-a_lambda.*Wn).*((Senn_plus_lambda_3*Lp/Dp).*cosh(Wn/Lp) + sinh(Wn/Lp)));
J_n_lambda_4 = -k1_lambda.*Lp.*a_lambda.*exp(-a_lambda.*Wn) + (k1_lambda./((Senn_plus_lambda_4.*Lp/Dp).*sinh(Wn/Lp)+cosh(Wn/Lp))).*((Senn_plus_lambda_4.*Lp/Dp)+ a_lambda.*Lp - exp(-a_lambda.*Wn).*((Senn_plus_lambda_4*Lp/Dp).*cosh(Wn/Lp) + sinh(Wn/Lp)));

J_p_lambda = k2_lambda .* a_lambda .* Ln - (k2_lambda./((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln))).*((Sb*Ln/Dn).*(cosh(Wp/Ln) - exp(-a_lambda.*Wp)) + sinh(Wp/Ln) + a_lambda.*Ln.*exp(-a_lambda*Wp)) ;
J_p_lambda_1 = k2_lambda_1 .* a_lambda .* Ln - (k2_lambda_1./((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln))).*((Sb*Ln/Dn).*(cosh(Wp/Ln) - exp(-a_lambda.*Wp)) + sinh(Wp/Ln) + a_lambda.*Ln.*exp(-a_lambda*Wp)) ;
J_p_lambda_2 = k2_lambda_2 .* a_lambda .* Ln - (k2_lambda_2./((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln))).*((Sb*Ln/Dn).*(cosh(Wp/Ln) - exp(-a_lambda.*Wp)) + sinh(Wp/Ln) + a_lambda.*Ln.*exp(-a_lambda*Wp)) ;
J_p_lambda_3 = k2_lambda_3 .* a_lambda .* Ln - (k2_lambda_3./((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln))).*((Sb*Ln/Dn).*(cosh(Wp/Ln) - exp(-a_lambda.*Wp)) + sinh(Wp/Ln) + a_lambda.*Ln.*exp(-a_lambda*Wp)) ;
J_p_lambda_4 = k2_lambda_4 .* a_lambda .* Ln - (k2_lambda_4./((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln))).*((Sb*Ln/Dn).*(cosh(Wp/Ln) - exp(-a_lambda.*Wp)) + sinh(Wp/Ln) + a_lambda.*Ln.*exp(-a_lambda*Wp)) ;

J_Wd_lambda = e .* f_lambda.*(1-R_lambda).*exp(-a_lambda*Wn).*(1-exp(-a_lambda*Wd)); %current of depletion zone 
J_tot_lambda = Jn0_plus_x1_lambda + J_Wa_x3_lambda +J_p_lambda +J_n_lambda +J_Wd_lambda;
J_tot_lambda_1 = Jn0_plus_x1_lambda_1+ J_Wa_x3_lambda_1 +J_p_lambda_1 +J_n_lambda_1 +J_Wd_lambda;
J_tot_lambda_2 =  Jn0_plus_x1_lambda_2 + J_Wa_x3_lambda_2 +J_p_lambda_2 +J_n_lambda_2 +J_Wd_lambda;
J_tot_lambda_3 =  Jn0_plus_x1_lambda_3 + J_Wa_x3_lambda_3 +J_p_lambda_3 +J_n_lambda_3 +J_Wd_lambda;
J_tot_lambda_4 =  Jn0_plus_x1_lambda_4 + J_Wa_x3_lambda_4 +J_p_lambda_4 +J_n_lambda_4 +J_Wd_lambda;
 figure;
 plot(lambda, J_tot_lambda,'r-', lambda, J_tot_lambda_1,'y-', lambda, J_tot_lambda_2,'g-', lambda, J_tot_lambda_3,'k-', lambda, J_tot_lambda_4,'b-');
 title('Total photocurrent density VS wavelength')
 legend('Wn_plus = 0.01e-6','Wn_plus = 0.02e-6','Wn_plus = 0.03e-6','Wn_plus = 0.04e-6','Wn_plus = 0.05e-6')
 ylabel('Photocurrent density A/m2');
 xlabel('wavelength (µm)');

 figure;

 plot(lambda, J_p_lambda,'b-', lambda, J_n_lambda, 'g-',lambda, J_Wd_lambda, 'r-',lambda, Jn0_plus_x1_lambda,'y-', lambda,J_Wa_x3_lambda,'k-' );
 title('Contribution to the photocurrent density')
 legend('contribution p region','contribution n region','contribution depletion zone', 'contribution over doped zone', 'contribution Wa');
  ylabel('Photocurrent density A/m2');
 xlabel('wavelength (µm)');
%Photocurent density en fonction de Wn
Senn_plus_Wn = (Nd/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %ne depend pas de Wn donc blc
Senn_plus_Wn_1 = (Nd_1/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %ne depend pas de Wn donc blc
Senn_plus_Wn_2 = (Nd_2/Nd_plus) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %ne depend pas de Wn donc blc

Sen_plus_n_Wn = (Nd_plus/Nd)*(Dp/Lp)*coth(Wn_range/Ln); %depend de Wn donc matrice
Sen_plus_n_Wn_1 = (Nd_plus/Nd_1)*(Dp/Lp)*coth(Wn_range/Ln); %depend de Wn donc matrice
Sen_plus_n_Wn_2 = (Nd_plus/Nd_2)*(Dp/Lp)*coth(Wn_range/Ln); %depend de Wn donc matrice

Fh_1_Wn = 1./(1+(Senn_plus_Wn./Sen_plus_n_Wn).*(Nd_plus/Nd));
Fh_1_Wn_1 = 1./(1+(Senn_plus_Wn_1./Sen_plus_n_Wn_1).*(Nd_plus/Nd_1));
Fh_1_Wn_2 = 1./(1+(Senn_plus_Wn_2./Sen_plus_n_Wn_2).*(Nd_plus/Nd_2));

Fh_1_corrected_Wn= Fh_1_Wn .*sech(Wn_range./Lp); %matrice
Fh_1_corrected_Wn_1= Fh_1_Wn_1 .*sech(Wn_range./Lp); %matrice
Fh_1_corrected_Wn_2= Fh_1_Wn_2 .*sech(Wn_range./Lp); %matrice

Tau_lambda_fixe = f(lambda_fixe)*((1-R(lambda_fixe))*e*Lp_plus*a(lambda_fixe))/(-1+(a(lambda_fixe)*Lp_plus)^2); 
Jn0_plus_x1_Wn = Fh_1_corrected_Wn.*(Tau_lambda_fixe * ((-a(lambda_fixe)*Lp_plus).*exp(-a(lambda_fixe)*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a(lambda_fixe)*Lp_plus - exp(-a(lambda_fixe)*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))));
Jn0_plus_x1_Wn_1 = Fh_1_corrected_Wn_1.*(Tau_lambda_fixe * ((-a(lambda_fixe)*Lp_plus).*exp(-a(lambda_fixe)*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a(lambda_fixe)*Lp_plus - exp(-a(lambda_fixe)*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))));
Jn0_plus_x1_Wn_2 = Fh_1_corrected_Wn_2.*(Tau_lambda_fixe * ((-a(lambda_fixe)*Lp_plus).*exp(-a(lambda_fixe)*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a(lambda_fixe)*Lp_plus - exp(-a(lambda_fixe)*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))));

J_Wa_x3_Wn = e*f(lambda_fixe)*(1-R(lambda_fixe))*exp(-a(lambda_fixe)*Wn_plus)*(1-exp(-a(lambda_fixe)*vpa(Wa)))*sech(Wn_range/Lp);

k1 = e * f(lambda_fixe) * ((1-R(lambda_fixe)) * a(lambda_fixe)* Lp/((a(lambda_fixe)*Lp)^2 -1 )); 
k2 = (e*f(lambda_fixe)*(1-R(lambda_fixe))*a(lambda_fixe)*Ln/(((a(lambda_fixe)*Ln)^2) -1))*exp(-a(lambda_fixe)*(Wn+Wn_plus+Wa+Wd));
J_n_Wn = -k1.*Lp.*a(lambda_fixe).*exp(-a(lambda_fixe).*Wn_range) + (k1./((Senn_plus_Wn.*Lp/Dp).*sinh(Wn_range/Lp)+cosh(Wn_range/Lp))).*((Senn_plus_Wn.*Lp/Dp)+ a(lambda_fixe).*Lp - exp(-a(lambda_fixe).*Wn_range).*((Senn_plus_Wn*Lp/Dp).*cosh(Wn_range/Lp) + sinh(Wn_range/Lp))); % photocurrent generated by p layer calculated at x3
J_p_Wn = k2 * a(lambda_fixe) * Ln - (k2/((Sb*Ln/Dn)*sinh(Wp/Ln) + cosh(Wp/Ln)))*((Sb*Ln/Dn)*(cosh(Wp/Ln) - exp(-a(lambda_fixe)*Wp)) + sinh(Wp/Ln) + a(lambda_fixe)*Ln*exp(-a(lambda_fixe)*Wp)) ;

J_Wd = e * f(lambda_fixe)*(1-R(lambda_fixe))*exp(-a(lambda_fixe)*Wn_range)*(1-exp(-a(lambda_fixe)*Wd));
J_tot = J_p_Wn + J_n_Wn  + Jn0_plus_x1_Wn + J_Wa_x3_Wn + J_Wd ;
J_tot_1 = J_p_Wn + J_n_Wn  + Jn0_plus_x1_Wn_1 + J_Wa_x3_Wn + J_Wd ;
J_tot_2 = J_p_Wn + J_n_Wn  + Jn0_plus_x1_Wn_2 + J_Wa_x3_Wn + J_Wd ;


figure;
hold on;
plot(Wn_range, J_n_Wn,'mx');

plot(Wn_range,J_p_Wn,'co');

plot(Wn_range, Jn0_plus_x1_Wn,'bs');

plot(Wn_range,J_Wa_x3_Wn,'gd');

plot(Wn_range, J_Wd,'kd');
legend('n contribution', 'p contribution','n+ contribution','Wa contribution', 'Wd contribution');
title('contribution to photocurrent VS Wn thickness');

xlabel('Wn (m)');
ylabel('Photocurrent A/m2');

figure;
plot(Wn_range, J_tot,'m-',Wn_range,J_tot_1,'c-', Wn_range,J_tot_2,'r-');
xlabel('Wn (m)');
 ylabel ('Photocurrent A/m2');
 title('Total photocurrent density');
 
 %Internal quantum efficiency VS lambda (we don(t take in account
 %reflectance and transmission across the cell
 P_in = 1000; % W/m2 according to AM = 1.5

nb_photons_impinging = 10^(-6)*P_in * lambda /(h * c); %nb photons *s-1*m-2 whose energy is h*c/lambda
nb_carriers_generated = J_tot_lambda  * 6.241e18 ;% np charges elementaires s-1 m-2
IQE_lambda = nb_carriers_generated./nb_photons_impinging;
figure;
plot(lambda,IQE_lambda);
title('Internal quantum efficiency VS wavelength');
xlabel('lambda (µm)');
ylabel('IQE');
%on fait varier Nd_plus pour afficher IQE en fct de lambda
Senn_plus_lambda_11 = (Nd/Nd_plus_1) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %effect of n+ on the n carriers 
Senn_plus_lambda_22= (Nd/Nd_plus_2) * (Dp_plus / Lp_plus) * ((Sp * Lp_plus / Dp_plus )+ tanh(Wn_plus / Lp_plus))/(1+ (Sp*Lp_plus/Dp_plus)* tanh(Wn_plus/Lp_plus)); %effect of n+ on the n carriers 

Sen_plus_n_lambda_11 = (Nd_plus_1/Nd)*(Dp/Lp)*coth(Wn/Ln); %effect of the n region on the n+ carriers 
Sen_plus_n_lambda_22 = (Nd_plus_2/Nd)*(Dp/Lp)*coth(Wn/Ln); %effect of the n region on the n+ carriers 

Fh_1_lambda_11 = 1/(1+(Senn_plus_lambda_11./Sen_plus_n_lambda_11)*(Nd_plus_1/Nd));
Fh_1_lambda_22 = 1/(1+(Senn_plus_lambda_22./Sen_plus_n_lambda_22)*(Nd_plus_2/Nd));

Fh_1_corrected_lambda_11= Fh_1_lambda_11 *sech(Wn/Lp); 
Fh_1_corrected_lambda_22= Fh_1_lambda_22 *sech(Wn/Lp); 

Jn0_plus_x1_lambda_11 = Fh_1_corrected_lambda_11*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))); 
Jn0_plus_x1_lambda_22 = Fh_1_corrected_lambda_22*Tau .* ((-a_lambda*Lp_plus).*exp(-a_lambda*Wn_plus) + (((Sp*Lp_plus/Dp_plus) + a_lambda*Lp_plus - exp(-a_lambda*Wn_plus)*((Sp*Lp_plus*cosh(Wn_plus/Lp_plus)/Dp_plus) + sinh(Wn_plus/Lp_plus)))/((Sp*Lp_plus*sinh(Wn_plus/Lp_plus)/Dp_plus) + cosh(Wn_plus/Lp_plus)))); 

J_tot_lambda_11 = Jn0_plus_x1_lambda_11 + J_Wa_x3_lambda +J_p_lambda +J_n_lambda +J_Wd_lambda;
J_tot_lambda_22=  Jn0_plus_x1_lambda_22 + J_Wa_x3_lambda +J_p_lambda +J_n_lambda +J_Wd_lambda;

nb_carriers_generated_11 = J_tot_lambda_11  * 6.241e18 ;% np charges elementaires s-1 m-2
nb_carriers_generated_22 = J_tot_lambda_22  * 6.241e18 ;% np charges elementaires s-1 m-2

IQE_lambda_11 = nb_carriers_generated_11./nb_photons_impinging;
IQE_lambda_22 = nb_carriers_generated_22./nb_photons_impinging;


plot(lambda,IQE_lambda, lambda, IQE_lambda_11, lambda,IQE_lambda_22);

%on fait varier Wn+ pour afficher IQE en fct de lambda

nb_carriers_generated_111 = J_tot_lambda_1*  6.241e18;
nb_carriers_generated_222 = J_tot_lambda_2*  6.241e18;
nb_carriers_generated_333 = J_tot_lambda_3 * 6.241e18;
nb_carriers_generated_444 = J_tot_lambda_4  * 6.241e18;

IQE_lambda_111 = nb_carriers_generated_111./nb_photons_impinging;
IQE_lambda_222 = nb_carriers_generated_222./nb_photons_impinging;
IQE_lambda_333 = nb_carriers_generated_333./nb_photons_impinging;
IQE_lambda_444 = nb_carriers_generated_444./nb_photons_impinging;

figure;
plot(lambda, IQE_lambda,lambda,IQE_lambda_111,lambda,IQE_lambda_222,lambda,IQE_lambda_333,lambda,IQE_lambda_444);










