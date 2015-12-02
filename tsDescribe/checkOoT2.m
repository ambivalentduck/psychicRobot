clc
clear all

N=10000;

OoT2=1+exprnd(1,N,4);
T=OoT2.^-.5;
gam=sum(T,2).^-2;

figure(1)
clf
ecdf(gam)
gamfit(gam)