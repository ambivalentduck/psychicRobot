clc
clear all
close all

A=5;

N=10000;

r=exprnd(A,N,2).*(sign(rand(N,2)-.5));

hist([sum(r,2),r(:,1),r(:,2)],20)

[counts,bins]=hist([sum(r,2),r(:,1),r(:,2)],40);
figure(2)
plot(abs(bins),log(counts),'.')