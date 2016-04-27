clc
clear all

M=1000;
lambda=1;

n=3.7;

summed=mod(n,1)*exprnd(lambda,M,1)+gamrnd(floor(n),lambda,M,1);

figure(1)
clf
hold on

ecdf(summed,'bounds','on')
[f,x]=ecdf(summed);
plot(x,gamcdf(x,n,lambda),'r')