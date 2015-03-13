function [ofinterest,Pnormecdf,Enormecdf]=boltzmannsubjectplots(t,x,v,a)

global fitmecounts fitmebins

%% Set constants
numbins=40;
lpctile=0; %This needs to be fixed instead by throwing out "stops"
upctile=99; %Throw away the craziest outliers because of Farnaz's coding quirks

sWidth=2;
sHeight=2;

%% Compute Pt, Tt
m=150*0.453592*.0597; %Dempster says arm is 5.97% of tbm, assume 150 lbs. Obviously effective mass is different and subject specific.
Pt=m*dot(v',a')';
Tt=.5*m*dot(v',v');

%% Bin and count Pt levels, plot raw and fit
subplot(sWidth,sHeight,1)
hold on

aPt=abs(Pt);
upper=min(prctile(aPt,upctile),10); %Power above ten watts is thoroughly implausible. Let's not fit bugs in Farnaz's code.
lower=prctile(aPt,lpctile);

edges=linspace(lower,upper,65)';
counts=zeros(64,1);
for k=2:length(edges)
    counts(k-1)=sum((aPt>=edges(k-1))&(aPt<edges(k)));
end
counts=counts/sum(counts);

centers=(edges(2:end)+edges(1:end-1))/2;
bar(centers,counts)

fitmecounts=counts;
fitmebins=edges;

W=fminbnd(@fitexpcdf,0,10);

ofinterest.W=W;

Pnormecdf=[counts centers/W];

cdf=expcdf(fitmebins,W);
fit=cdf(2:end)-cdf(1:end-1);
ofinterest.power_rmse=sqrt(mean((fit-counts).^2));
plot(centers,fit,'r','linewidth',3)

title(['RMS Error = ',num2str(ofinterest.power_rmse,3)])
xlabel('Power, Watts')
ylabel('Fraction of Time Points')
legend('Data',['Exponential Fit, W=',num2str(W,2)])

%% Check symmetry by reflecting across zero
edges=linspace(lower,upper,65)';
centers=(edges(2:end)+edges(1:end-1))/2;
pcounts=zeros(64,1);
ncounts=zeros(64,1);
nPt=-Pt;
for k=2:length(edges)
    pcounts(k-1)=sum((Pt>=edges(k-1))&(Pt<edges(k)));
    ncounts(k-1)=sum((nPt>=edges(k-1))&(nPt<edges(k)));
end
pcounts=pcounts/sum(pcounts);
ncounts=ncounts/sum(ncounts);

ofinterest.power_symmetry=sqrt(mean((pcounts-ncounts).^2));

subplot(sWidth,sHeight,2)
plot(centers,pcounts,'b',centers,ncounts,'r')
title(['RMS Error = ',num2str(ofinterest.power_symmetry,3)])
xlabel('Power, Watts')
ylabel('Fraction of Time Points')
legend('Positive Power','Negative Power')


%% Bin and count Tt levels, plot raw and fit
subplot(sWidth,sHeight,3)
hold on

upper=prctile(Tt,upctile);
lower=prctile(Tt,lpctile);

edges=linspace(lower,upper,65)';
counts=zeros(64,1);
for k=2:length(edges)
    counts(k-1)=sum((Tt>=edges(k-1))&(Tt<edges(k)));
end
counts=counts/sum(counts);

centers=(edges(2:end)+edges(1:end-1))/2;
bar(centers,counts)

fitmecounts=counts;
fitmebins=edges;
nT=fmincon(@fitgamcdf,[1.5; 1],-eye(2),-[1; 0.01]);

cdf=gamcdf(fitmebins,nT(1),nT(2));
fit=cdf(2:end)-cdf(1:end-1);
ofinterest.kinetic_rmse=sqrt(mean((fit-counts).^2));
plot(centers,fit,'r','linewidth',3)

title(['RMS Error = ',num2str(ofinterest.kinetic_rmse,3)])
xlabel('Kinetic Energy, Joules')
ylabel('Fraction of Time Points')
legend('Data',['Gamma Fit, n=',num2str(nT(1),2),', T=',num2str(nT(2),2)])

%% Bin and count speeds, plot raw and fit
subplot(sWidth,sHeight,4)
hold on

speed=vecmag(v);

upper=prctile(speed,upctile);
lower=prctile(speed,lpctile);

ofinterest.maxT=upper; %Because we have some CRAZY outliers which are completely incredible

edges=linspace(lower,upper,65)';
counts=zeros(64,1);
for k=2:length(edges)
    counts(k-1)=sum((speed>=edges(k-1))&(speed<edges(k)));
end
counts=counts/sum(counts);

centers=(edges(2:end)+edges(1:end-1))/2;
bar(centers,counts)

fitmecounts=counts;
fitmebins=.5*m*edges.^2;

nT=fmincon(@fitgamcdf,[1.5; 1],-eye(2),-[1; 0.01]);
cdf=gamcdf(fitmebins,nT(1),nT(2));
fit=cdf(2:end)-cdf(1:end-1);
ofinterest.n=nT(1);
ofinterest.T=nT(2);
ofinterest.speed_rmse=sqrt(mean((fit-counts).^2));

Enormecdf=[counts centers/nT(2)];

cdf=gamcdf(fitmebins,nT(1),nT(2));
fit=cdf(2:end)-cdf(1:end-1);

plot(centers,fit,'r','linewidth',3)

title(['RMS Error = ',num2str(ofinterest.speed_rmse,3)])
xlabel('Speed, m/s')
ylabel('Fraction of Time Points')
legend('Data','Gamma Fit')


end

