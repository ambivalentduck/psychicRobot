clc
clear all

N=10000;

OoT2=.5+exprnd(1,N,4);
gam=sum(OoT2,2);

figure(1)
clf
[x,y,z]=fitShiftedGam(gam,1)

figure(2)
clf
[F,x]=ecdf(gam);
plot(x,smooth(gradient(F)./gradient(x)))