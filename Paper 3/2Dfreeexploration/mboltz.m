function [p,a]=mboltz(x,a)

if nargin<2
    a=mean(x)/(2*sqrt(2/pi));
end
x2=x.^2;
p=.01*sqrt(2/pi)*x2.*exp(-x2./(2*a.^2))./(a.^3); %Wikipedia was off by a factor of 100?