clc
clear all

%profile on

addpath ../../DeOpt/

load ../Data/curlkick/curlkick1Y.mat

doPlots=1;
alpha=5;

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))

if doPlots
    figure(1002)
    clf
    hold on
end

alllumps=0;

for ff=6 %1:length(f)
    T=f(ff);
    yoff=ff-1;
    
    y=trials(T).y;
    t=trials(T).ty;
    inds=trials(T).i0:trials(T).if+20;
    [~,peaks]=findpeaks(-vecmag(y(inds,3:4)));
    inds=inds(1):inds(peaks(end)-5);
    t=t-t(inds(1));
    
    [lumps,de(ff).resid]=findBetaLumps(t',y(:,3:4),inds,alpha);
    gtlumps=lumps;
    de(ff).lumps=lumps;
    plot(t',vecmag(y(:,3:4))+yoff,'k')
    plot(t',vecmag(de(ff).resid)+yoff,'k-.')
    re(ff).y=zeros(length(trials(T).y),2);
    
    for k=1:length(lumps)
        tau=(t'-lumps(k).C)/lumps(k).S+.5;
        tau=max(min(tau,1),0);
        kappa=betapdf(tau,alpha,alpha)/lumps(k).S;
        re(ff).y=re(ff).y+kappa*lumps(k).L;
        if doPlots
            plot(t',vecmag(kappa*lumps(k).L)+yoff,'b')
        end
    end
    if doPlots
        plot(t,vecmag(re(ff).y)+yoff,'c')
    end
    
    [lumps,re(ff).resid]=findBetaLumps(t',re(ff).y,inds);
    %lumps=deoptSubunits(lumps,t(inds)',re(ff).y(inds,:));
    de(ff).y=zeros(length(trials(T).y),2);
    if doPlots
        plot(t',vecmag(re(ff).resid)+yoff,'r-.')
    end
    for k=1:length(lumps)
        tau=(t'-lumps(k).C)/lumps(k).S+.5;
        tau=max(min(tau,1),0);
        kappa=betapdf(tau,alpha,alpha)/lumps(k).S;
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

%profile viewer

%% Plot for Paper/Konrad

figure(1006)
clf
set(gcf,'color','w')

subplot(2,2,1)
hold on
ul=[0 2];
plot(ul,ul,'m')
plot(cdiffP(:,2),cdiffP(:,1),'.','markersize',.75)
axis equal
set(gca,'ytick',[])
xlim(ul)
ylim(ul)
xlabel('C, seconds')
set(gca,'xaxisLocation','top')
set(gca,'ycolor','w')

subplot(2,2,2)
hold on
ul=[0 1];
plot(ul,ul,'m')
plot(sdiffP(:,2),sdiffP(:,1),'.','markersize',.75)
axis equal
set(gca,'ytick',[])
xlim(ul)
ylim(ul)
xlabel('S, seconds')
set(gca,'xaxisLocation','top')
set(gca,'ycolor','w')

subplot(2,2,3)
hold on
ul=[-.25 .25];
plot(ul,ul,'m')
plot(lxdiffP(:,2),lxdiffP(:,1),'.','markersize',.75)
axis equal
set(gca,'ytick',[])
set(gca,'xtick',[-.2 .2])
xlim(ul)
ylim(ul)
xlabel('L_x, m/s')
set(gca,'ycolor','w')

subplot(2,2,4)
hold on
ul=[-.25 .25];
plot(ul,ul,'m')
plot(lydiffP(:,2),lydiffP(:,1),'.','markersize',.75)
axis equal
set(gca,'ytick',[])
set(gca,'xtick',[-.2 .2])
xlim(ul)
ylim(ul)
xlabel('L_y, m/s')
set(gca,'ycolor','w')

subplot(2,2,1)
set(gca,'position',[0 .5 .4 .4])
subplot(2,2,2)
set(gca,'position',[.35 .5 .4 .4])
subplot(2,2,3)
set(gca,'position',[0 .1 .4 .4])
subplot(2,2,4)
set(gca,'position',[.35 .1 .4 .4])