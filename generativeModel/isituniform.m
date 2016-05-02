clc
clear all

M=1000;

Nmin=2;
Nave=3.88;
theta=Nave-Nmin;
ns=Nmin+round(exprnd(theta,M,1));

[f,x]=ecdf(ns);
figure(1)
clf
plot(x,log(1-f))

maxes=zeros(M,1);
for k=1:M
    %vals(k).x=rand(1,ns(k));
    vals(k).x=exprnd(1,1,ns(k));
    vals(k).x=1*vals(k).x/sum(vals(k).x);
    maxes(k)=max(vals(k).x);
end

vcat=[vals.x];
[f,x]=ecdf(vcat);
figure(2)
clf
subplot(2,1,1)
plot(x,log(1-f),'k.')
subplot(2,1,2)
[f,x]=ecdf(maxes);
plot(x,f)