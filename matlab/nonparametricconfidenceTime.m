clc
clear all

for k=1 %:4
    load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])

    clustermeS=zeros(length(trials)-1,1);
    clustermeE=clustermeS;
    for t=1:length(trials)-1
        clustermeS(t)=trials(t+1).x(1,1);
        clustermeE(t)=trials(t+1).x(end,1);
    end
    [cats,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
    means=sort(means);

    %For each direction/length/disturbance (2*2*5), plot an example from subject S
    %UNDER that example, plot the t-statistic: (Y-mean)/(sd(baseline)*sqrt(1+1/Nbaseline)
    % 95% confidence t<=.2

    nclean=2;
    clean=0*dcats;
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
    end
    clean=(clean==0)';

    figure(k)
    clf
    hold on

    bins=0:.006:1; %Should give about the same number of total bins, but we can adjust this empirically

    for ST=1:3
        for EN=1:3
            if ST==EN
                continue
            end
            f=find((starts==ST)&(ends==EN)&clean);
            [ST EN length(f)]
            
            figure(ST*10+EN)
            clf
            hold on
            for kf=1:length(f)
                speed=vecmag(trials(f(kf)).v);
                onset=find(speed>.05,1,'first');
                start=max(onset-35,1);
                plot(trials(f(kf)).t-trials(f(kf)).t(start),speed)
            end
        end
    end
end

