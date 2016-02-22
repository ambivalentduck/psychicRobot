clc
clear all

profile on

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
    
    [lumps(ff).fit,de(ff).resid]=findLumps(t',y(:,3:4),inds);
    
    resid=y(:,3:4);
    
    for k=1 %:length(lumps(ff).fit)
        C=lumps(ff).fit(k).C;
        S=lumps(ff).fit(k).S;
        L=lumps(ff).fit(k).L;
        
        [~,ipeak]=min(abs(t-C));
        [~,ilower]=min(abs(t-C+S/2));
        [~,iupper]=min(abs(t-C-S/2));
        Vpeak=resid(ipeak,:);
        lumps(ff).raw(k).v=S*resid(ilower:iupper,:)/mean(resid(ilower:iupper,:));
        lumps(ff).raw(k).t=(t(ilower:iupper)-t(ipeak))/S+.5;
        plot(lumps(ff).raw(k).t,lumps(ff).raw(k).v,'b')
        %plot(lumps(ff).raw(k).t,vecmag(resid(ilower:iupper,:)),'r.','markersize',.75)
        
        tau=(t'-C)/S+.5;
        tau=max(min(tau,1),0);
        kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/S;
        resid=resid-kappa*L;
    end   
end

rawcat=[lumps.raw];
tcat=[rawcat.t]';
vcat=vertcat(rawcat.v);

minme = @(n) sum((2*vcat-betapdf(tcat,n,n)).^2);

n=fminunc(minme,1.5)
plot(0:.01:1,.5*betapdf(0:.01:1,n,n),'k-','linewidth',3)

