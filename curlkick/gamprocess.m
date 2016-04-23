clc
clear all
close all

load lumps.mat
colors=linspecer(8);
figure(1)
msize=10;

hold on
for k=1:8
    V=vertcat(lumps(:,k).L2).*vertcat(lumps(:,k).Tn2);
    [fV,xV]=ecdf(V);
    p=polyfit(xV(1:end-1),log(1-fV(1:end-1)),1);
    
    n=vertcat(lumps(:,k).n);
    [fn,xn]=ecdf(n);
    plot(xV,log(1-fV),'.','color',colors(k,:),'markersize',msize)
    
end
