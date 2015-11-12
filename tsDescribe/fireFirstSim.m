clc
clear all

N=1000;
CANDIDATES=10000;
A=9;

chosen=zeros(N,1);
for k=1:N
    ISI=exprnd(A,CANDIDATES,1);
    chosen(k)=min(ISI);
end

figure(1)
clf
hold on
mu=expfit(chosen);
[f,x,flo,fhi]=ecdf(chosen);
h=fill([x(~isnan(flo)); wrev(x(~isnan(fhi)))],[flo(~isnan(flo)); wrev(fhi(~isnan(fhi)))],'b');
set(h,'edgealpha',0,'facealpha',.5)
plot(x,f,'b')
plot(x,expcdf(x,mu),'r')

mu-A/CANDIDATES