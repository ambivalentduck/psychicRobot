clc
clear all
close all

load lumps.mat
colors=linspecer(8);
figure(1)
msize=10;

subplot(2,2,1)
hold on
for k=1:8
    subs(k).resid=vertcat(lumps(:,k).resid);
    subs(k).resid=subs(k).resid/mean(subs(k).resid);
    [f,x]=ecdf(subs(k).resid);
    plot(x,log(1-f),'.','color',colors(k,:),'markersize',msize)
    subs(k).residf=f;
    subs(k).residx=x;
end
R2(1)=linfitR2(vertcat(subs.residx),vertcat(subs.residf));
title('Residuals')

subplot(2,2,2)
hold on
for k=1:8
    subs(k).L2Tn2=vertcat(lumps(:,k).L2).*vertcat(lumps(:,k).Tn2);
    subs(k).L2Tn2=subs(k).L2Tn2/mean(subs(k).L2Tn2);
    [f,x]=ecdf(subs(k).L2Tn2);
    plot(x,log(1-f),'.','color',colors(k,:),'markersize',msize)
    subs(k).L2Tn2f=f;
    subs(k).L2Tn2x=x;
end
R2(2)=linfitR2(vertcat(subs.L2Tn2x),vertcat(subs.L2Tn2f));
title('Peak Kinetic Energy')

subplot(2,2,4)
hold on
nDCs=zeros(8,3);
for k=1:8
    subs(k).n=vertcat(lumps(:,k).n);
    [f,x]=ecdf(subs(k).n);
    inds=find(~(0==(1-f)));
    lrf=log(1-f(inds));
    x=x(inds);
    p=polyfit(x,lrf,1);
    nDCs(k,1:2)=[k p(1)];
    subs(k).n=subs(k).n/mean(subs(k).n);
    [f,x]=ecdf(subs(k).n);
    plot(x,log(1-f),'.','color',colors(k,:),'markersize',msize)
    subs(k).nf=f;
    subs(k).nx=x;
end
R2(3)=linfitR2(vertcat(subs.nx),vertcat(subs.nf));
title('Submotion Count')

subplot(2,2,3)
hold on
for k=1:8
    subs(k).dc=abs(horzcat(lumps(:,k).dC)').^2;
    nDCs(k,3)=mean(subs(k).dc);
    subs(k).dc=subs(k).dc/mean(subs(k).dc);
    [f,x]=ecdf(subs(k).dc);
    
    plot(x,log(1-f),'.','color',colors(k,:),'markersize',msize)
    %p=gamfit(subs(k).dc)
    %plot(x,gamcdf(x,p(1),p(2)),'-','color',colors(k,:))
    
    subs(k).dcf=f;
    subs(k).dcx=x;
end
%R2(4)=linfitR2(vertcat(subs.dcx),vertcat(subs.dcf))
title('Peak-to-Peak Duration')
nDCs


for k=1:4
    subplot(2,2,k)
    %xlim([0 4])
    %ylim([-6 0])
    text(4,0,['R^2 = ',num2str(R2(k),2)],'horizontalalignment','right','verticalalignment','top')
    %set(gca,'xticklabels',{'0','2u
end
suplabel('log(1-Cumulative Probability)','y')
suplabel('Mean Normalized Value','x')

