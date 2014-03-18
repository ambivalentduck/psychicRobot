clc
clear all

figure(5)
clf
hold on

colors=[1 0 0; 0 1 0; 0 0 1; 1 1 1]; %RGBK

S=3;
ANECDOTES=ones(2*2*5,1);

for k=1 %:4
    load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])

    %For each direction/length/disturbance (2*2*5), plot an example from subject S
    %UNDER that example, plot the t-statistic: (Y-mean)/(sd(baseline)*sqrt(1+1/Nbaseline)
    % 95% confidence t<=.2

    figure(1)
    clf
    hold on
    clustermeS=zeros(length(trials)-1,1);
    clustermeE=clustermeS;
    for t=1:length(trials)-1
        plot(trials(t+1).x([1 end],1),trials(t+1).x([1 end],2),'.')
        clustermeS(t)=trials(t+1).x(1,1);
        clustermeE(t)=trials(t+1).x(end,1);
    end
    [cats,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
    centers=sort(means);
    plot(means,.5*[1 1 1],'rx')
    axis equal

    direction=sign(starts-ends);
    distance=abs(starts-ends);

    reachcat=10*starts+ends;
    urc=unique(reachcat);
    urc=urc(2:end);

    for U=1:length(urc)
        st=floor(urc(U)/10);
        en=mod(urc(U),10);
        edges=[-inf min(0,centers(en)-centers(st)):.0015:max(0,centers(en)-centers(st)) inf]; % 1.5mm bins
        for kk=1:length(baselineCatme(U,:))
            if ~isempty(baselineCatme(U,kk).x)
                baselineCatme(U,kk).x(1,:)=baselineCatme(U,kk).x(1,:)-baselineCatme(U,kk).x(1,1);
                baselineCatme(U,kk).y(1,:)=baselineCatme(U,kk).y(1,:)-baselineCatme(U,kk).y(1,1);
            end
        end

        X=[baselineCatme(U,:).x];
        Y=[baselineCatme(U,:).y];
        baselines(U).means=zeros(length(edges)-1,1);
        baselines(U).stds=zeros(length(edges)-1,1);
        baselines(U).counts=zeros(length(edges)-1,1);
        baselines(U).edges=edges;
        for kk=1:length(edges)-1
            inds=find((X(1,:)>=edges(kk))&(X(1,:)<edges(kk+1)));
            baselines(U).means(kk)=mean(Y(1,inds));
            baselines(U).stds(kk)=std(Y(1,inds));
            baselines(U).counts(kk)=length(inds);
        end
    end
    figure(6)
    clf
    hold on
    N=20;
    for U=1:length(urc)
        plot(baselines(U).edges(2:end-1),baselines(U).means(2:end)+U/N)
        plot(baselines(U).edges(2:end-1),baselines(U).means(2:end)-baselines(U).stds(2:end)+U/N,'-.')
        plot(baselines(U).edges(2:end-1),baselines(U).means(2:end)+baselines(U).stds(2:end)+U/N,'-.')
    end
    axis equal

    %Organize into columns: short, short, long, long
    %Organize into rows: disturbances 1:5

    outsidethelines=0;
    figure(678)
    clf
    hold on

    for dir=[-1 1]
        for dist=[1 2]
            %This actually contains TWO reachcats if dist==1
            for dcat=1:5
                F=find((direction==dir)&(distance==dist)&(dcats'==dcat));
                for kk=1:length(F)
                    c=F(kk);
                    rc=find(urc==(10*starts(c)+ends(c)));
                    xoff=(dir+3*dist)/3; %maps to 2 4 5 7
                    yoff=dcat/20;

                    Y=trials(c).y(:,2);
                    Y=Y-Y(1);
                    X=trials(c).y(:,1);
                    X=X-X(1);

                    if kk==1 %ANECDOTES(
                        figure(5)
                        plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+yoff)
                        plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)-1.96*baselines(rc).stds(2:end)+yoff,'-.')
                        plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+1.96*baselines(rc).stds(2:end)+yoff,'-.')
                        plot(X+xoff,Y+yoff,'r')
                    end

                    tstat=zeros(length(baselines(rc).edges)-1,1);

                    for kkk=1:length(baselines(rc).edges)-1
                        inds=find((X>=baselines(rc).edges(kkk))&(X<baselines(rc).edges(kkk+1)));
                        if isempty(inds)
                            tstat(kkk)=nan;
                        else
                            tstat(kkk)=abs(mean(Y(inds))-baselines(rc).means(kkk))/(baselines(rc).stds(kkk)); %*sqrt(1+1/baselines(rc).counts(kkk)));
                        end
                    end
                    myo2=mean(baselines(rc).means);
                    figure(5)
                    for kkk=2:length(baselines(rc).edges)-1
                        if tstat(kkk)>=1.96
                            plot(baselines(rc).edges(kkk)+xoff,myo2-.01-.01*(kk/length(F))+yoff,'.','color',colors(k,:),'markersize',.001)
                        end
                    end

                    if sum(tstat)>0
                        outsidethelines=outsidethelines+1
                        figure(678)
                        xoff=0;
                        yoff=outsidethelines/10;
                        plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+yoff)
                        plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)-1.96*baselines(rc).stds(2:end)+yoff,'.')
                        plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+1.96*baselines(rc).stds(2:end)+yoff,'.')
                        plot(X+xoff,Y+yoff,'r-o')
                        plot(trials(c).x(:,1)-trials(c).x(1,1),trials(c).x(:,2)-trials(c).x(1,2)+yoff,'g-x')
                        tst=tstat(2:end);
                        tst(tst>.2)=.25;
                        plot(baselines(rc).edges(2:end-1)+xoff,tst/10-.03+yoff,'k')
                    end
                end

            end
        end
    end
end
axis equal