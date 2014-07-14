clc
clear all

d=0:.001:.1;
o=.02;
syms b real
a=double(solve(1/(1+exp(-b*(-.01)))-.05))

val=-1;
s=1+(val-1)*(1./(1+exp(-a*(d-o))));

figure(1)
clf
plot(d,s)