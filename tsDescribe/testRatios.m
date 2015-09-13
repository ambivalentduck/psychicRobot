clc
clear all

N=100000;
A=5;

L2=A+gamrnd(3,A,N,1);
T2=exprnd(A,N,1);

rat=L2.*T2;

figure(1)
clf
hold on
ecdf(rat)
mu=expfit(rat)
fitline=0:.1:50;
plot(fitline,expcdf(fitline,mu),'r')