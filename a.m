function a = absorption(lambda) 
if (0.88 <= lambda ) 
    a = 0;
 
elseif (0.8 <= lambda ) && (lambda < 0.88)
    a = 10^(-37.5*lambda +34)*100; 
elseif (0.8 > lambda)
    a = 10^(-3.3*lambda + 6.64)*100;
    
end


end