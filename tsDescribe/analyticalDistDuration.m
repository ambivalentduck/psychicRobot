clc
clear all

syms a n theta lambda positive
syms x real

gpdf=1/(gamma(n)*theta^n)*x^(n-1)*exp(-x/theta);

apdf=symsum(1/lambda*exp(-n/lambda)*gpdf,n,2,inf)

pretty(apdf)