clc
clear all
close all

load lumps.mat
colors=linspecer(8);
figure(1)
msize=10;

means=zeros(8,1);

hold on
for k=1:8
    first=zeros(45,1)-1;
    last=first;
    for kk=1:45
        if length(lumps(kk,k).dC)<2
            continue
        end
        first(kk)=lumps(kk,k).dC(1);
        last(kk)=lumps(kk,k).dC(end);
    end
    first=first(first>=0);
    last=last(last>=0);
    
    subplot(1,2,1)
    hold on
    dC=[lumps(:,k).dC];
    means(k)=mean(dC.^2);
    [f,x]=ecdf(dC.^2./mean(dC.^2));
    plot(x,log(1-f),'color',colors(k,:))
    
    subplot(1,2,2)
    hold on
    
    [ff,fx]=ecdf(first.^2);
    [lf,lx]=ecdf(last.^2);
    plot(fx,log(1-ff),'-','color',colors(k,:),'markersize',msize)
    plot(lx,log(1-lf),'.','color',colors(k,:),'markersize',msize)
    continue
    p=gamfit(first)
    plot(fx,gamcdf(fx,p(1),p(2)),'-','color',colors(k,:))
end
xlim([0 1])
legend('first','last')
xlabel('Mean-Normalized Inter-Submotion Interval Squared')
subplot(1,2,1)
ylabel('Cumulative Probability')
xlabel('Inter-Submotion Peak Interval Squared')

means
mean(means)
