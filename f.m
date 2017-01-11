function f = flux(lambda) 
if (0.24 <= lambda ) && (lambda <= 0.48) 
    f = 100*100*1.5*(19.7*lambda -4.7)*10^15;
 
elseif (0.48 < lambda ) && (lambda <= 1.08)
    f = 100*100*1.5*(-2.5*lambda + 5.7)*10^15;
end 


end 