function [shift,n,T]=fitTs(SUB)

load(['../Data/Data_pulse/pulse',num2str(SUB),'.mat'])

%% Categorize by start/end pair
clustermeS=zeros(length(trials)-1,1); %#ok<NODEF>
clustermeE=clustermeS;
for t=1:length(trials)-1
    clustermeS(t)=trials(t+1).x(1,1);
    clustermeE(t)=trials(t+1).x(end,1);
end
[trash,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
means=sort(means);

starts=zeros(length(trials),1);
ends=starts;
for t=1:length(trials)
    [trash,starts(t)]=min(abs(trials(t).x(1,1)-means));
    [trash,ends(t)]=min(abs(trials(t).x(end,1)-means));
    trialInfo(t).startcat=starts(t);
    trialInfo(t).endcat=ends(t);
    trialInfo(t).x0=means(starts(t));
    trialInfo(t).xf=means(ends(t));
end
starts(1)=-inf; %never use this reach...
trialInfo(1).startcat=-inf;

reachcat=10*starts+ends;
urc=unique(reachcat);
urc=urc(2:end);

%% Use undisturbed examples only
dcats=[trials.disturbcat];

nclean=4;
clean=0*dcats;
for kk=1:nclean
    clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
end
clean=(clean==0);

for t=1:length(trials)
    trialInfo(t).clean=clean(t);
    trialInfo(t).reachcat=find(reachcat(t)==urc);
end

%% For each trial, fit a 5th order poly to velocity

tesses=-1*ones(size(trials));

for t=2:length(trials)
    if trialInfo(t).clean
        tesses(t)=getTsMetric(trials(t).x,trials(t).t,[means(starts(t)) .5],[means(ends(t)) .5]);
    end
end

figure(SUB)
clf
subplot(1,2,1)
hold on
tessesfixed=abs(tesses(tesses~=-1));
[shift,n,T]=fitShiftedGam(tessesfixed);
x=tessesfixed.^-2;
x=x(x<35); %Not censoring, just visibility for plotting

nbins=25;
edges=linspace(min(x),max(x),nbins+1)';
counts=zeros(nbins,1);
for k=2:length(edges)
    counts(k-1)=sum((x>=edges(k-1))&(x<=edges(k)));
end
counts=counts/sum(counts);

centers=(edges(2:end)+edges(1:end-1))/2;
bar(centers,counts)
gcx=gamcdf(edges-shift,n,T);
plot(centers,gcx(2:end)-gcx(1:end-1),'r')
ylabel('Normalized Frequency')
xlabel('t_s^{-2}')
title('Histogram of t_s^{-2}')
subplot(1,2,2)
hold on
ecdf(x,'bounds','on')
sortX=sort(x);
plot(sortX,gamcdf(sortX-shift,n,T),'r')
title(['U=',num2str(shift),' n=',num2str(n),' T=',num2str(T)])

