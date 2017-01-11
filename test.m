clear all


for i=1:100
    x(i)=i;
    f(i)=sin(1/x(i));
end

figure; plot(x,f);
