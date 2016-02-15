clc
clear all

profile on

addpath ../../DeOpt/

load ../Data/curlkick/curlkick1Y.mat

doPlots=1;

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))

if doPlots
    figure(1002)
    clf
    hold on
end

alllumps=0;

for ff=1:length(f)
    T=f(ff);
    yoff=ff-1;
    
    y=trials(T).y;
    t=trials(T).ty;
    inds=trials(T).i0:trials(T).if+20;
    [~,peaks]=findpeaks(-vecmag(y(inds,3:4)));
    inds=inds(1):inds(peaks(end)-5);
    t=t-t(inds(1));
    
    [lumps,de(ff).resid]=findLumps(t',y(:,3:4),inds);
    gtlumps=lumps;
    de(ff).lumps=lumps;
    plot(t',vecmag(y(:,3:4))+yoff,'k')
    plot(t',vecmag(de(ff).resid)+yoff,'k-.')
    re(ff).y=zeros(length(trials(T).y),2);
    
    for k=1:length(lumps)
        tau=(t'-lumps(k).C)/lumps(k).S+.5;
        tau=max(min(tau,1),0);
        kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/lumps(k).S;
        re(ff).y=re(ff).y+kappa*lumps(k).L;
        if doPlots
            plot(t',vecmag(kappa*lumps(k).L)+yoff,'b')
        end
    end
    if doPlots
        plot(t,vecmag(re(ff).y)+yoff,'c')
    end
    
    [lumps,re(ff).resid]=findLumps(t',re(ff).y,inds);
    %lumps=deoptSubunits(lumps,t(inds)',re(ff).y(inds,:));
    de(ff).y=zeros(length(trials(T).y),2);
    if doPlots
        plot(t',vecmag(re(ff).resid)+yoff,'r-.')
    end
    for k=1:length(lumps)
        tau=(t'-lumps(k).C)/lumps(k).S+.5;
        tau=max(min(tau,1),0);
        kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/lumps(k).S;
        de(ff).y=de(ff).y+kappa*lumps(k).L;
        if doPlots
            plot(t',vecmag(kappa*lumps(k).L)+yoff,'r--')
            drawnow
        end
    end
    if doPlots
        plot(t,vecmag(de(ff).y)+yoff,'m')
    end
    
    for k=1:length(lumps)
        alllumps=alllumps+1;
        C=lumps(k).C;
        [~,ic]=min(abs([de(ff).lumps.C]-C));
        cdiff(alllumps)=lumps(k).C-de(ff).lumps(ic).C;
        lxdiff(alllumps)=lumps(k).L(1)-de(ff).lumps(ic).L(1);
        lydiff(alllumps)=lumps(k).L(2)-de(ff).lumps(ic).L(2);
        sdiff(alllumps)=lumps(k).S-de(ff).lumps(ic).S;
        cdiffP(alllumps,1:2)=[lumps(k).C,de(ff).lumps(ic).C];
        lxdiffP(alllumps,1:2)=[lumps(k).L(1),de(ff).lumps(ic).L(1)];
        lydiffP(alllumps,1:2)=[lumps(k).L(2),de(ff).lumps(ic).L(2)];
        sdiffP(alllumps,1:2)=[lumps(k).S,de(ff).lumps(ic).S];
    end
end

figure(1005)
clf
subplot(2,2,1)
ecdf(cdiff)
xlabel('C error, seconds')
subplot(2,2,2)
ecdf(sdiff)
xlabel('S error, seconds')
subplot(2,2,3)
ecdf(lxdiff)
xlabel('L_x error, meters per second')

subplot(2,2,4)
ecdf(lydiff)
xlabel('L_y error, meters per second')

profile viewer

figure(1006)
clf
subplot(2,2,1)
plot(cdiffP(:,1),cdiffP(:,2),'.')
xlabel('C error, seconds')
subplot(2,2,2)
plot(sdiffP(:,1),sdiffP(:,2),'.')
xlabel('S error, seconds')
subplot(2,2,3)
plot(lxdiffP(:,1),lxdiffP(:,2),'.')
xlabel('L_x error, meters per second')
subplot(2,2,4)
plot(lydiffP(:,1),lydiffP(:,2),'.')
xlabel('L_y error, meters per second')
