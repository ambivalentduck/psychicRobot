clc
clear all
close all

%Set ranges
%For each range, show the bar plots with figure window titled
%    set(gcf,'name',range)
%Grab the mean and confidence interval from c=multcompare with output off
%Plot mean and confidence interval against range
%call that finalfig5

NR=12;

for k=1:NR
    ranges{k}=(0:4)+(k-1)*5;
    f=figure(k);
    hold on
end

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'U.mat'])

    %In theory, you need to sort by reachcat, which is in reachinfo

    %For now ignore that loop, just sort by subject by baseline, extracted
    %baseline, disturbed hand, extracted from disturbance full and partial,
    %so 5 data per subject.

    clear cleanStruct
    f=find([trialInfo.clean]);
    for t=2:length(f)
        lt=length(trials(f(t)).t);
        lty=length(undisturbed(f(t)).ty);
        inds=(1+lt-lty):lt;
        cleanStruct(t).R2BX=getR2(trials(f(t)).x(inds,:));
        cleanStruct(t).R2BY=getR2(undisturbed(f(t)).y);
    end
    subject(k,1).R2BX=vertcat(cleanStruct.R2BX);
    subject(k,1).R2BY=vertcat(cleanStruct.R2BY);

    clear triStruct

    for R=1:length(ranges)
        RANGE=ranges{R};
        for t=1:length(trials)
            if isempty(trialInfo(t).forceinds)
                thishappens=1
                continue
            end
            inds=trialInfo(t).forceinds(1)+RANGE;
            onset=find(vecmag(trials(t).v)>.05,1,'first');
            start=max(onset-35,1);
            indsy=max(1,trialInfo(t).forceinds(1)-start)+RANGE;

            triStruct(t).R2X=getR2(trials(t).x(inds,:));
            triStruct(t).R2Y=getR2(trials(t).y(indsy,:));
        end
        subject(k).R2X=vertcat(triStruct.R2X);
        subject(k).R2Y=vertcat(triStruct.R2Y);

        CBX=1*ones(size(subject(k).R2BX));
        CBY=2*ones(size(subject(k).R2BY));
        CX=3*ones(size(subject(k).R2X));
        CY=4*ones(size(subject(k).R2Y));

        Xlist=[CBX;CBY;CX;CY];
        Ylist=[subject(k).R2BX;subject(k).R2BY;subject(k).R2X;subject(k).R2Y];

        xoff=.3/2*((k-1)-1.5);
        figure(R)
        plot(Xlist(:,1)+xoff+.05*(rand(size(Xlist(:,1)))-.5),Ylist(:,1),'g.','markersize',1)
        boxplot(Ylist(:,1),Xlist,'orientation','vertical','notch','on','positions',(1:4)+xoff,'widths',.1,'symbol','r.','labels',{'','','',''})

        lists(k,R).x=Xlist(:,1);
        lists(k,R).s=Xlist(:,1)*0+k;
        lists(k,R).y=Ylist(:,1);
    end
end

for k=1:NR
    figure(k)
    set(gcf,'name',[num2str(5*ranges{k}(1)),'-',num2str(5*ranges{k}(end)),' ms'])
    ylabel('Max Perpendicular Deviation, cm')
    set(gca,'xtick',[1:4])
    set(gca,'xticklabels',{'Hand','Extracted','Hand','Extracted'})
    xlim([.7 4.3])
    ylim([0 17])
    ys=ylim;
    plot(2.5*[1 1],ys,'m-')
    ylim(ys)

    text(1.5,-1,'Disturbed Reaches','horizontalalignment','center')
    text(3.5,-1,'Undisturbed Reaches','horizontalalignment','center')
end

figure(8) %Have to know figure of interest off-hand
set(0,'defaulttextinterpreter','none')
set(gcf,'color',[1 1 1])
set(gcf,'position',[625   250   440   490])
laprint(gcf,'figures/fig4raw','scalefonts','off','asonscreen','on')

figure(NR+1)

return

fullX=vertcat(lists.x);
fullS=vertcat(lists.s);
fullY=vertcat(lists.y);
[p,table,stats]=kruskalwallis(fullY,fullX)
multcompare(stats)
f=find(fullX==1);
%[p,table,stats]=kruskalwallis(fullY(f),fullS(f))
%multcompare(stats)

f2=find(fullX==2);
mean(fullY(f2)-fullY(f))
std(fullY(f2)-fullY(f))

f3=find(fullX==3);
f4=find(fullX==4);

mean(fullY(f3)-fullY(f4))
std(fullY(f3)-fullY(f4))

