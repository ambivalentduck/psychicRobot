clc
clear all
close all

global binedges bincounts nT

sublets={'JL','JT','LP','MF','YM'};

nbins=64;
lpctile=0; %This needs to be fixed instead by throwing out "stops"
upctile=99; %Throw away the craziest outliers because of Farnaz's coding quirks

sH=4;
sW=1;

for k=4
    fname=['free_exp_',sublets{k},'.mat'];
    load(fname)
    figure(k)
    clf

    f=find((vecmag(v)>.2)|(vecmag(a)>1));
    t=t(f);
    x=x(f,:);
    v=v(f,:);
    a=a(f,:);

    m=180*0.453592*.0597; %Dempster says arm is 5.97% of tbm, assume 180 lbs. Obviously effective mass is different and subject specific.
    Pt=m*dot(v',a')';
    Tt=.5*m*dot(v',v');

    subplot(sH,sW,1)
    hold on

    aPt=abs(Pt);
    upper=min(prctile(aPt,upctile),10); %Power above ten watts is thoroughly implausible. Let's not fit bugs or hiccups.
    lower=prctile(aPt,lpctile);
    edges=linspace(lower,upper,nbins+1)';
    counts=zeros(nbins,1);
    for k=2:length(edges)
        counts(k-1)=sum((aPt>=edges(k-1))&(aPt<edges(k)));
    end
    counts=counts/sum(counts);
    centers=(edges(2:end)+edges(1:end-1))/2;
    bar(centers,counts)

    binedges=edges;
    bincounts=counts;
    lambda=fmincon(@fitexpcdf,1,-1,0)
    [cost,fit]=fitexpcdf(lambda);
    plot(centers,fit)

    subplot(sH,sW,2)
    hold on

    upper=min(prctile(Tt,upctile),1); %Power above ten watts is thoroughly implausible. Let's not fit bugs or hiccups.
    lower=prctile(Tt,lpctile);
    edges=linspace(lower,upper,nbins+1)';
    counts=zeros(nbins,1);
    for k=2:length(edges)
        counts(k-1)=sum((Tt>=edges(k-1))&(Tt<edges(k)));
    end
    counts=counts/sum(counts);
    centers=(edges(2:end)+edges(1:end-1))/2;
    bar(centers,counts)

    binedges=edges;
    bincounts=counts;
    nT=fmincon(@fitgamcdf,[1,1],-eye(2),-[1; .01])
    [cost,fit]=fitgamcdf(nT);
    plot(centers,fit,'r')

    subplot(sH,sW,3)
    hold on
    meanX=mean(x);
    tolx=.02;
    f=find(vecmag([x(:,1)-meanX(1) x(:,2)-meanX(2)])<tolx);
    Tf=Tt(f);

    counts=zeros(nbins,1);
    for k=2:length(edges)
        counts(k-1)=sum((Tf>=edges(k-1))&(Tf<edges(k)));
    end
    counts=counts/sum(counts);
    centers=(edges(2:end)+edges(1:end-1))/2;
    bar(centers,counts)
    binedges=edges;
    bincounts=counts;
    offset=fminunc(@fitgamoffset,0)
    [cost,fit]=fitgamoffset(offset);
    plot(centers,fit,'r')

    xgrid=min(x(:,1)):tolx:max(x(:,1));
    ygrid=min(x(:,2)):tolx:max(x(:,2));
    offset=zeros(length(xgrid),length(ygrid));
    warning off all
    for xn=1:length(xgrid)
        for yn=1:length(ygrid)
            f=find(vecmag([x(:,1)-xgrid(xn) x(:,2)-ygrid(yn)])<tolx);
            Tf=Tt(f);
            counts=zeros(nbins,1);
            for k=2:length(edges)
                counts(k-1)=sum((Tf>=edges(k-1))&(Tf<edges(k)));
            end
            counts=counts/sum(counts);

            binedges=edges;
            bincounts=counts;
            offset(xn,yn)=fminunc(@fitgamoffset,0);
        end
    end
    warning on all
end

%Pick a pixel, say...mean x +/- .5 cm
figure(75)
clf
hold on
plot3(x(:,1),x(:,2),Tt,'k.','markersize',1)

[mxgrid,mygrid]=meshgrid(xgrid,ygrid);
mesh(mxgrid,mygrid,offset')



