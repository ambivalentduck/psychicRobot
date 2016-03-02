clc
clear all

global fit_t fit_y lumpsvec

load ../Data/curlkick/curlkick1Y.mat

clean_init=1;

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))

alllumps=0;
hulls=zeros(length(f),1);
cleanhulls=zeros(length(f),1);

for ff=1:length(f)
    T=f(ff);
    yoff=ff-1;
    
    %% Load and format trial
    y=trials(T).y;
    t=trials(T).ty;
    inds=trials(T).i0:trials(T).if+20;
    [~,peaks]=findpeaks(-vecmag(y(inds,3:4)));
    inds=inds(1):inds(peaks(end)-5);
    t=t-t(inds(1));
    
    
    [lumps,de(ff).resid]=findLumps(t',y(:,3:4),inds);
    [~,hulls(ff)]=convhull(de(ff).resid);
    
    %% Optional polish on S, L
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
    
    %% Reconstitute the reach AND perform metrics starting from the end and building up
    de(ff).lumps=lumps;
    re(ff).y=zeros(length(trials(T).y),2);
    re(ff).y=re(ff).y+de(ff).resid
    for k=length(lumps):-1:1
        y=lumpy(t,lumps(k));
        re(ff).y=re(ff).y+y;
        lump=findLumps(t,re(ff).y,inds,1);
        
        
        alllumps=alllumps+1;

        cdiff(alllumps)=lump.C-lumps(k).C;
        ldiff(alllumps)=norm(lump.L-lumps(k).L);
        sdiff(alllumps)=lump.S-lumps(k).S;
        cdiffP(alllumps,1:2)=[lump.C,lumps(k).C];
        lxdiffP(alllumps,1:2)=[lump.L(1),lumps(k).L(1)];
        lydiffP(alllumps,1:2)=[lump.L(2),lumps(k).L(2)];
        sdiffP(alllumps,1:2)=[lump.S,lumps(k).S];
    end
    
    if clean_init
        [~,cleanhulls(ff)]=convhull(fit_y-re(ff).y(inds,:));
    end
    
    
end

figure(1005)
clf
subplot(2,2,1)
ecdf(cdiff)
xlabel('C error, seconds')
subplot(2,2,2)
ecdf(sdiff)
hold on
%[f,x]=ecdf(sdiffPS(:,1)-sdiffPS(:,2));
%plot(x,f,'r')
xlabel('S error, seconds')
subplot(2,2,3)
ecdf(ldiff)

%% Plot for Paper/Konrad

figure(1006)
clf
set(gcf,'color','w')
msize=10;

subplot(2,2,1)
hold on
ul=[0 2];
plot(ul,ul,'m')
plot(cdiffP(:,2),cdiffP(:,1),'.','markersize',msize)
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