clc
clear all
close all

load lumps.mat

colors=.2+.6*rand(8,3);

figure(1)
hold on

for k=1:8
    subs(k).resid=vertcat(lumps(:,k).resid);
    subs(k).resid=subs(k).resid/mean(subs(k).resid);
    [f,x]=ecdf(subs(k).resid);
    plot(x,log(1-f),'.','color',colors(k,:))
end

resid=vertcat(subs.resid);
[f,x]=ecdf(resid);
plot(x,log(1-f),'-')