clc
clear all

global fit_t fit_y lumpsvec

load ../Data/curlkick/curlkick1Y.mat

doPlots=0;
clean_init=1;

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))

if doPlots
    figure(1002)
    clf
    hold on
end

alllumps=0;
coslumpN=0;
hulls=zeros(length(f),1);
nlumps=zeros(length(f),1);

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
    nlumps(ff)=length(lumps);

    [~,hulls(ff)]=convhull(de(ff).resid);
    if clean_init
        fit_y=y(inds,3:4);
        fit_t=t(inds);
        lumpsvec=lumps2vec(lumps);
        warning off all
        Sopt=fminunc(@justSobj,[lumps.S]',optimset('TolFun',1e-12));
        warning on all
        [justSobj([lumps.S]'),justSobj(Sopt),norm([lumps.S]'-Sopt)]
        lumps=vec2lumps(adjustLumpS(lumpsvec,Sopt));
    end
    
    de(ff).lumps=lumps;
    re(ff).y=zeros(length(trials(T).y),2);
    re(ff).ycos=zeros(length(trials(T).y),2);
    
    if doPlots
        figure(1002)
        plot(t',vecmag(y(:,3:4))+yoff,'k')
        plot(t',vecmag(de(ff).resid)+yoff,'k-.')
    end
    
    for k=1:length(lumps)
        y=lumpy(t,lumps(k));
        yc=cosLumpy(t,lumps(k));
        re(ff).y=re(ff).y+y;
        re(ff).ycos=re(ff).ycos+yc;
        if doPlots
            plot(t',vecmag(y)+yoff,'b')
        end
    end
    if doPlots
        plot(t,vecmag(re(ff).y)+yoff,'c')
        plot(t,vecmag(re(ff).ycos)+yoff,'c-.')
    end
    if clean_init
        [~,cleanhulls(ff)]=convhull(fit_y-re(ff).y(inds,:));
    end
    
    [lumps,re(ff).resid]=findLumps(t',re(ff).y,inds);
    [lumpscos,re(ff).resid]=findLumps(t',re(ff).ycos,inds,5);
    if doPlots
        figure(1002)
    end
    %lumpSopt=deoptSubunits(lumps2vec(lumps),t(inds)',re(ff).y(inds,:));
    lumpsvec=lumps2vec(lumps);
    fit_y=re(ff).y;
    fit_t=t;
    warning off all
    Sopt=fminunc(@justSobj,[lumps.S]',optimset('TolFun',1e-16));
    warning on all
    [justSobj([lumps.S]'),justSobj(Sopt),norm([lumps.S]'-Sopt)]
    lumpSopt=vec2lumps(adjustLumpS(lumpsvec,Sopt));
    
    de(ff).y=zeros(length(trials(T).y),2);
    deS(ff).y=zeros(length(trials(T).y),2);
    if doPlots
        plot(t',vecmag(re(ff).resid)+yoff,'r-.')
    end
    
    for k=1:length(lumps)
        rey=lumpy(t,lumps(k));
        yopt=lumpy(t,lumpSopt(k));
        de(ff).y=de(ff).y+rey;
        deS(ff).y=deS(ff).y+yopt;
        if doPlots
            plot(t',vecmag(rey)+yoff,'r--')
            plot(t',vecmag(yopt)+yoff,'g--')
            drawnow
        end
    end
    if doPlots
        plot(t,vecmag(de(ff).y)+yoff,'m')
        plot(t,vecmag(deS(ff).y)+yoff,'y')
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
        sdiffPS(alllumps,1:2)=[lumpSopt(k).S,de(ff).lumps(ic).S];
    end
    for k=1:length(lumpscos)
        coslumpN=coslumpN+1;
        C=lumpscos(k).C;
        [~,ic]=min(abs([de(ff).lumps.C]-C));
        cdiffCos(coslumpN)=lumpscos(k).C-de(ff).lumps(ic).C;
        lxdiffCos(coslumpN)=lumpscos(k).L(1)-de(ff).lumps(ic).L(1);
        lydiffCos(coslumpN)=lumpscos(k).L(2)-de(ff).lumps(ic).L(2);
        sdiffCos(coslumpN)=lumpscos(k).S-de(ff).lumps(ic).S;
        cdiffPCos(coslumpN,1:2)=[lumpscos(k).C,de(ff).lumps(ic).C];
        lxdiffPCos(coslumpN,1:2)=[lumpscos(k).L(1),de(ff).lumps(ic).L(1)];
        lydiffPCos(coslumpN,1:2)=[lumpscos(k).L(2),de(ff).lumps(ic).L(2)];
        sdiffPCos(coslumpN,1:2)=[lumpscos(k).S,de(ff).lumps(ic).S];
    end
end

figure(1005)
clf
subplot(2,2,1)
hold on
ecdf(cdiff)
[f,x]=ecdf(cdiffCos);
plot(x,f,'g')
xlabel('C error, seconds')
subplot(2,2,2)
ecdf(sdiff)
hold on
[f,x]=ecdf(sdiffPS(:,1)-sdiffPS(:,2));
plot(x,f,'r')
[f,x]=ecdf(sdiffCos);
plot(x,f,'g')
xlabel('S error, seconds')
subplot(2,2,3)
ecdf(lxdiff)
hold on
[f,x]=ecdf(lxdiffCos);
plot(x,f,'g')
xlabel('L_x error, meters per second')

subplot(2,2,4)
ecdf(lydiff)
hold on
[f,x]=ecdf(lydiffCos);
plot(x,f,'g')
xlabel('L_y error, meters per second')

%% Plot for Paper/Konrad

figure(1006)
clf
set(gcf,'color','w')
msize=2;

subplot(2,2,1)
hold on
ul=[0 2];
plot(ul,ul,'m')
plot(cdiffP(:,2),cdiffP(:,1),'.','markersize',msize)
plot(cdiffPCos(:,2),cdiffPCos(:,1),'g.','markersize',msize)
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
plot(sdiffP(:,2),sdiffP(:,1),'b.','markersize',msize)
plot(sdiffPCos(:,2),sdiffPCos(:,1),'g.','markersize',msize)
plot(sdiffPS(:,2),sdiffPS(:,1),'r.','markersize',msize)
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
plot(lxdiffP(:,2),lxdiffP(:,1),'.','markersize',msize)
plot(lxdiffPCos(:,2),lxdiffPCos(:,1),'g.','markersize',msize)
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
plot(lydiffP(:,2),lydiffP(:,1),'.','markersize',msize)
plot(lydiffPCos(:,2),lydiffPCos(:,1),'g.','markersize',msize)
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

figure(2006)
clf
hold on
%[f,x]=ecdf(hulls);
%plot(x,f,'b')
ecdf(cleanhulls,'bounds','on')
[f,x]=ecdf(cleanhulls);
mu=expfit(cleanhulls);
plot(x,expcdf(x,mu),'r')