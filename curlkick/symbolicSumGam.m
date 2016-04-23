clc
clear all

syms n x theta lambda a positive

gampd=1/(gamma(n)*theta^n)*x^(n-1)*exp(-x/theta);

symsum(gampd*1/lambda*exp(-n/lambda),n,1,inf)