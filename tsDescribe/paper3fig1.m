clc
clear all

t=0:.001:1;
tc=.5;
ts=1;
ta=(t-tc)/ts+.5;
K=(30*ta.^2-60*ta.^3+30*ta.^4)/ts;

figure(1)
clf
hold on
ecdf(K)
vRange=0:.001:1.875;
Pv=1-sqrt(1-2*sqrt(30)*sqrt(vRange)/15);

Erange=0:.001:(.5*1.875^2);
PE=1-sqrt(1-2/15*sqrt(30)*sqrt(sqrt(2)*sqrt(y)));
plot(out,Pout,'r.')

%Make random walk