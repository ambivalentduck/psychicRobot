function reachStruct=processPulse(SUB)

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
    reachStruct(t).x=trials(t).x;
    reachStruct(t).v=trials(t).v;
    reachStruct(t).t=trials(t).t;
    reachStruct(t).startcat=starts(t);
    reachStruct(t).endcat=ends(t);
    reachStruct(t).x0=[means(starts(t)) .5];
    reachStruct(t).xf=[means(ends(t)) .5];
end
starts(1)=-inf; %never use this reach...
reachStruct(1).startcat=-inf;

reachcat=10*starts+ends;
urc=unique(reachcat);
urc=urc(2:end);

%% Use undisturbed examples only
dcats=[trials.disturbcat];

nclean=2;
clean=0*dcats;
for kk=1:nclean
    clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
end
clean=(clean==0);
clean(1)=0;

for t=1:length(trials)
    reachStruct(t).clean=clean(t);
    reachStruct(t).reachcat=find(reachcat(t)==urc);
end

reachStruct=reachStruct([reachStruct.clean]);

if nargout<1
    figure(57)
    clf
    hold on
    for t=1:length(reachStruct)
        plot(reachStruct(t).x(:,1),reachStruct(t).v(:,1),'linewidth',.1)
    end
    
    fitTsArray(reachStruct)
end
