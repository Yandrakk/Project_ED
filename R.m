function R = reflectance(lambda) 
%n = sqrt(5.372514 + (5.466742*lambda^2)/(lambda^2 - 0.4431307^2) + (0.02429960*lambda^2)/(lambda^2 - 0.8746453^2) + (1.957522*lambda^2)/(lambda^2 -36.9166^2));
%R = abs((n-1)/(n+1))^2;
if (lambda <= 0.25)
    R = 3*lambda -0.15;
elseif(lambda > 0.25)
    R  = -0.3017*lambda + 0.675;
end
 
end
