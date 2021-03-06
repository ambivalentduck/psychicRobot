function ofinterest=boltzmannsubjectplots(t,x,v,a)

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
startW=mean(aPt(aPt>=edges(1))&(aPt<edges(end)));
expfit=fmincon(@bruteforceexp,1,-1,-1E-5);

ofinterest.W=expfit;

cdf=expcdf(fitmebins,expfit);
fit=cdf(2:end)-cdf(1:end-1);
ofinterest.power_rmse=sqrt(mean((fit-counts).^2));
plot(centers,fit,'r','linewidth',3)

title(['RMS Error = ',num2str(ofinterest.power_rmse,3)])
xlabel('Power, Watts')
ylabel('Fraction of Time Points')
legend('Data',['Exponential Fit, W=',num2str(expfit,2)])

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
gamfit=fmincon(@bruteforcegam,[1.5; .5] ,-eye(2),[-1;-0.01]);

cdf=gamcdf(fitmebins,gamfit(1),gamfit(2));
fit=cdf(2:end)-cdf(1:end-1);
ofinterest.kinetic_rmse=sqrt(mean((fit-counts).^2));
plot(centers,fit,'r','linewidth',3)

title(['RMS Error = ',num2str(ofinterest.kinetic_rmse,3)])
xlabel('Kinetic Energy, Joules')
ylabel('Fraction of Time Points')
legend('Data',['Gamma Fit, n=',num2str(gamfit(1),2),', T=',num2str(gamfit(2),2)])

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
gamfit=fmincon(@bruteforcegam,[1.5; .5] ,-eye(2),[-1;-0.01]);
ofinterest.n=gamfit(1);
ofinterest.T=gamfit(2);

cdf=gamcdf(fitmebins,gamfit(1),gamfit(2));
fit=cdf(2:end)-cdf(1:end-1);
ofinterest.speed_rmse=sqrt(mean((fit-counts).^2));
plot(centers,fit,'r','linewidth',3)

title(['RMS Error = ',num2str(ofinterest.speed_rmse,3)])
xlabel('Speed, m/s')
ylabel('Fraction of Time Points')
legend('Data',['Gamma Fit, n=',num2str(gamfit(1),2),', T=',num2str(gamfit(2),2)])


end

function cost=bruteforcegam(x)
global fitmecounts fitmebins

cdf=gamcdf(fitmebins,x(1),x(2));
fit=cdf(2:end)-cdf(1:end-1);
cost=sum((fitmecounts-fit).^2);

end

function cost=bruteforceexp(x)
global fitmecounts fitmebins

cdf=expcdf(fitmebins,x);
fit=cdf(2:end)-cdf(1:end-1);

cost=sum((fitmecounts-fit).^2);

end

