clc
clear all
close all

load lumps.mat
colors=linspecer(8);
figure(1)
msize=10;


hold on
for k=1:8
    subs(k).dc=horzcat(lumps(:,k).dC)';
    subs(k).dc=abs(subs(k).dc); %/mean(subs(k).dc);
    [f,x]=ecdf(subs(k).dc);
    if any(x<0)
        continue
    end
    plot(x,f,'.','color',colors(k,:),'markersize',msize)
    [shift,n,T]=fitShiftedGam(subs(k).dc);
    ns(k)=n;
    plot(x,gamcdf(x-shift,n,T),'-','color',colors(k,:))
    
    subs(k).dcf=f;
    subs(k).dcx=x;
end

n