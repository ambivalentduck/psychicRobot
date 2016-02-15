clc
clear all

load ../Data/curlkick/curlkick1Y.mat

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))

figure(1001)
clf
hold on

for ff=1:length(f)
    T=f(ff);
    
    y=trials(T).y;
    t=trials(T).ty;
    inds=trials(T).i0:trials(T).if;
    [~,peaks]=findpeaks(-vecmag(y(inds,:)));
    inds=inds(1):inds(peaks(end)-5);
    
    [lumps,resid]=findLumps(t',y(:,3:4),inds);
    
    colors=[1 0 0;
        0 1 0;
        0 0 1;
        .8 .8 0;
        0 1 1;
        1 0 1;
        rand(10,3)];
    
    if 0
    figure(ff)
    clf
    hold on
    plot(t(inds),vecmag(resid(inds,:)),'k')
    for k=1:length(lumps)
        tau=(t(inds)'-lumps(k).C)/lumps(k).S+.5;
        tau=max(min(tau,1),0);
        kappa=30*tau.^2-60*tau.^3+30*tau.^4;
        plot(t(inds)',vecmag(kappa*lumps(k).L),'color',colors(k,:))
    end
    end
    
    figure(1001)
    Y=cumtrapz(t(inds),resid(inds,:),1);
    plot(Y(:,1),Y(:,2)+ff/100,'k')
    axis equal
    [~,V(ff)] = convhull(resid(:,1),resid(:,2)); %Convex hull and its volume corresponding to an [X Y] set.
    [~,X(ff)] = convhull(1000*Y(:,1),1000*Y(:,2)); %Convex hull and its volume corresponding to an [X Y] set.
end

figure(1002)
clf
subplot(2,1,1)
hold on
[f,x]=ecdf(V);
ecdf(V,'bounds','on')
mu=gamfit(V);
plot(x,gamcdf(x,mu(1),mu(2)),'r')

subplot(2,1,2)
hold on
[f,x]=ecdf(V);
ecdf(V,'bounds','on')
mu=expfit(V);
plot(x,expcdf(x,mu),'r')