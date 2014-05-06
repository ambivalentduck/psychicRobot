clc
clear all
close all

%Set ranges
%For each range, show the bar plots with figure window titled
%    set(gcf,'name',range)
%Grab the mean and confidence interval from c=multcompare with output off
%Plot mean and confidence interval against range
%call that finalfig5

NR=21;

for k=1:NR
    ranges{k}=(0:4)+(k-1)*5;
end

green=[.1 .7 .3];

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'U.mat'])

    flag=0;

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
            if ~sum(trials(t).disturbcat==[1:4])
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


        if ~flag&(ranges{R}(3)*5>185)
            if k==1
                figure('name',[num2str(5*ranges{R}(1)),'-',num2str(5*ranges{R}(end)),' ms']);
                hold on
            end
            plot(Xlist(:,1)+xoff+.05*(rand(size(Xlist(:,1)))-.5),Ylist(:,1),'.','color',green,'markersize',2.5)
            boxplot(Ylist(:,1),Xlist,'orientation','vertical','notch','on','positions',(1:4)+xoff,'widths',.1,'symbol','r.','labels',{'','','',''})
            flag=1;
        end

        lists(k,R).x=Xlist(:,1);
        lists(k,R).s=Xlist(:,1)*0+k;
        lists(k,R).y=Ylist(:,1);
    end
end

h=findall(0,'type','figure');
for k=1:length(h)
    figure(h)
    %ylabel('Max Perpendicular Deviation, cm')
    plot(.64*[1 1],.1+[0 2],'k','linewidth',2)
    text(.61,1+.1,'2 cm','rotation',90,'horizontalalignment','center','verticalalignment','bottom')
    set(gca,'xtick',[1:4])
    set(gca,'xticklabel',{'Hand','Extracted','Hand','Extracted'})
    xlim([.6 4.4])
    plot([.7 4.3],[0 0],'k')
    ylim([0 14])
    %plot(2.5*[1 1],ys,'m-')
    %ylim(ys)
    text(1.5,-1,'Undisturbed Reaches','horizontalalignment','center')
    text(3.5,-1,'Disturbed Reaches','horizontalalignment','center')
    set(gca,'Box','off')
    set(gca,'YColor',[1 1 1])
    set(gca,'YTick',[])
    %set(gca,'XColor',[1 1 1])
end


for k=1:NR
    %for k=1:NR
    fullX=vertcat(lists(:,k).x);
    fullY=vertcat(lists(:,k).y);
    [p,table,stats]=anova1(fullY,fullX,'off');
    %[c{k} m{k}]=multcompare(stats,'display','off');
    mids(k,:)=stats.means;
    errors(k,:)=halfWidth(stats)';
    rangemids(k)=ranges{k}(1)*5;
end

set(0,'defaulttextinterpreter','none')
set(gcf,'color',[1 1 1])
set(gcf,'position',[625   250   440   490])
laprint(gcf,'figures/fig4raw','scalefonts','off','asonscreen','on')

red=[1 .3 .3];
green=[.1 .7 .3];

figure(NR+1)
clf
hold on

%fill([rangemids rangemids(end) rangemids(1)],[mids(:,2)+errors(:,2); 0; 0],'w','facecolor',green,'edgecolor',green)
%fill([rangemids rangemids(end) rangemids(1)],[mids(:,1)+errors(:,1); 0; 0],'w','facecolor',green,'edgecolor',red)
plot(rangemids,mids(:,1)+errors(:,1),'color',green)

for k=[3 4]
        if rem(k,2)~=0
            color=[0 0 0];
        else
            color=red;
        end


        fill([rangemids rangemids(end:-1:1)],[mids(:,k)-errors(:,k); mids(end:-1:1,k)+errors(end:-1:1,k)],'w','facecolor',color,'edgecolor','w')
        plot(rangemids,mids(:,k),'-','color',color)
end
xlabel('Time Post Onset of Disturbing Forces, ms')
ylabel('Max Perpendicular Distance, cm')




