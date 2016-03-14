clc
clear all

global fit_t fit_y lumpsvec

load ../Data/curlkick/curlkick1Y.mat

doPlots=0;
clean_init=0;

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))

alllumps=0;
coslumpN=0;

cdiff=[];
sdiff=[];
ldiff=[];
cdiffCos=[];
sdiffCos=[];
ldiffCos=[];

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

    [~,hulls(ff)]=convhull(de(ff).resid);

    de(ff).L2=sum(vertcat(lumps.L).^2,2);
    de(ff).Sn2=vertcat(lumps.S).^-2;
    de(ff).dC=diff(sort(vertcat(lumps.C)));
    de(ff).T=de(ff).L2.*de(ff).Sn2;
    de(ff).n=length(lumps);
    de(ff).lumps=lumps;
    re(ff).y=zeros(length(trials(T).y),2);
    re(ff).ycos=zeros(length(trials(T).y),2);
        
    for k=1:length(lumps)
        y=lumpy(t,lumps(k));
        re(ff).y=re(ff).y+y;
        yc=cosLumpy(t,lumps(k));
        re(ff).ycos=re(ff).ycos+yc;
    end
    
    [lumps,re(ff).resid]=findLumps(t',re(ff).y,inds);
    [lumpscos,re(ff).resid]=findLumps(t',re(ff).ycos,inds);
    re(ff).n=length(lumps);
    re(ff).ncos=length(lumpscos);
    
    for k=1:length(lumps)
        C=lumps(k).C;
        [~,ic]=min(abs([de(ff).lumps.C]-C));
        cdiff(end+1,:)=[lumps(k).C,de(ff).lumps(ic).C];
        sdiff(end+1,:)=[lumps(k).S,de(ff).lumps(ic).S];
        ldiff(end+1,:)=[lumps(k).L(1),de(ff).lumps(ic).L(1)];
        ldiff(end+1,:)=[lumps(k).L(1),de(ff).lumps(ic).L(1)];
    end
    
    for k=1:length(lumpscos)
        C=lumpscos(k).C;
        [~,ic]=min(abs([de(ff).lumps.C]-C));
        cdiffCos(end+1,:)=[lumpscos(k).C,de(ff).lumps(ic).C];
        sdiffCos(end+1,:)=[lumpscos(k).S,de(ff).lumps(ic).S];
        ldiffCos(end+1,:)=[lumpscos(k).L(1),de(ff).lumps(ic).L(1)];
        ldiffCos(end+1,:)=[lumpscos(k).L(1),de(ff).lumps(ic).L(1)];
    end
end

ndiff=[[de.n]',[re.n]'];
ndiffCos=[[de.n]',[re.ncos]'];

%% Plot for Paper/Konrad

figure(1)
clf
set(gcf,'color','w')
msize=3;

subplot(2,2,1)
hold on
ul=[0 2];
plot(ul,ul,'m')
plot(cdiff(:,2),cdiff(:,1),'k.','markersize',msize)
plot(cdiffCos(:,2),cdiffCos(:,1),'b.','markersize',msize)
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
plot(sdiff(:,2),sdiff(:,1),'k.','markersize',msize)
plot(sdiffCos(:,2),sdiffCos(:,1),'b.','markersize',msize)
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
plot(ldiff(:,2),ldiff(:,1),'k.','markersize',msize)
plot(ldiffCos(:,2),ldiffCos(:,1),'b.','markersize',msize)
axis equal
set(gca,'ytick',[])
set(gca,'xtick',[-.2 .2])
xlim(ul)
ylim(ul)
xlabel('L, m/s')
set(gca,'ycolor','w')

subplot(2,2,4)
hold on
ul=[2 16];
plot(ul,ul,'m')
plot(ndiff(:,2),ndiff(:,1),'k.','markersize',msize)
plot(ndiffCos(:,2),ndiffCos(:,1),'b.','markersize',msize)
axis equal
set(gca,'ytick',[])
set(gca,'xtick',[2 9 16])
xlim(ul)
ylim(ul)
xlabel('n')
set(gca,'ycolor','w')

subplot(2,2,1)
set(gca,'position',[0 .5 .4 .4])
subplot(2,2,2)
set(gca,'position',[.35 .5 .4 .4])
subplot(2,2,3)
set(gca,'position',[0 .1 .4 .4])
subplot(2,2,4)
set(gca,'position',[.35 .1 .4 .4])

%% Stats figure

figure(2)
clf
subplot(1,3,1)
hold on
n=ndiff(:,1);
[f,x]=ecdf(n);
plot(x,log(1-f),'k.')
ylabel('log(1-Cumulative Density)')
xlabel('n')

subplot(1,3,2)
hold on
%[f,x]=ecdf(hulls);
%plot(x,f,'b')
T=vertcat(de.T);
[f,x]=ecdf(T);
plot(x,log(1-f),'k.')
xlim([0 .1])
xlabel('L^2S^{-2}')

subplot(1,3,3)
hold on
%[f,x]=ecdf(hulls);
%plot(x,f,'b')
dC=vertcat(de.dC);
[f,x]=ecdf(dC);
plot(x,log(1-f),'k.')
xlim([0 .6])
xlabel('\Delta C')


