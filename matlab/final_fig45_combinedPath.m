clc
clear all
close all

%Set ranges
%For each range, show the bar plots with figure window titled
%    set(gcf,'name',range)
%Grab the mean and confidence interval from c=multcompare with output off
%Plot mean and confidence interval against range
%call that finalfig5

NR=10;
edges=linspace(0,100,NR+1); %sampling rate is 200 Hz = .005 s samples

black=[0 0 0];
red=[1 .3 .3];
green=[.1 .7 .3];

figure(1)
clf
hold on

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'U.mat'])

    flag=0;

    clear cleanStruct
    f=find([trialInfo.clean]);
    for t=2:length(f)
        lt=length(trials(f(t)).t);
        lty=length(undisturbed(f(t)).ty);
        inds=(1+lt-lty):lt; %Onset of movement
        cleanStruct(t).R2BX=getMaxPerp(trials(f(t)).x(inds,:));
        cleanStruct(t).R2BY=getMaxPerp(undisturbed(f(t)).y);
    end
    subject(k).R2BX=vertcat(cleanStruct.R2BX);
    subject(k).R2BY=vertcat(cleanStruct.R2BY);

    clear triStruct


    RANGE=edges(R):edges(R+1);
    for t=1:length(trials)
        if ~sum(trials(t).disturbcat==[1:4]) %Fast-ish skipping instead of indexing
            continue
        end
        inds=trialInfo(t).forceinds(1);
        onset=find(vecmag(trials(t).v)>.05,1,'first');
        start=max(onset-35,1);
        indsy=max(1,trialInfo(t).forceinds(1)-start);

        triStruct(t).R2X=getR2(trials(t).x(inds(1):end,:));
        triStruct(t).R2Y=getR2(trials(t).y(indsy(1);end,:));
    end

    for R=1:NR
        subject(k).R2X=vertcat(triStruct.R2X);
        subject(k).R2Y=vertcat(triStruct.R2Y);

        CBX=1*ones(size(subject(k).R2BX));
        CBY=4*ones(size(subject(k).R2BY));
        CX=2*ones(size(subject(k).R2X));
        CY=3*ones(size(subject(k).R2Y));

        %Xlist=[CBX;CBY;CX;CY];
        %Ylist=[subject(k).R2BX;subject(k).R2BY;subject(k).R2X;subject(k).R2Y];
        Xlist=[CBX;CX;CY];
        Ylist=[subject(k).R2BX;subject(k).R2X;subject(k).R2Y];


        DOTSIZE=4;
        LINEWIDTH=5;
        binCenter=5*(edges(R)+edges(R+1))/2;
        drawnWidth=5*.25*(edges(R+1)-edges(R));
        SCATTER=.2*drawnWidth;
        NLines=4;
        xoff=drawnWidth/2*(2/(NLines-1)*(k-1)-1);

        %         pctiles = prctile(subject(k).R2BX,[25;50;75]);
        %         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',green,'LineWidth',LINEWIDTH)
        %         pctiles = prctile(subject(k).R2X,[25;50;75]);
        %         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',black,'LineWidth',LINEWIDTH)
        %         pctiles = prctile(subject(k).R2Y,[25;50;75]);
        %         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',red,'LineWidth',LINEWIDTH)

        plot(zeros(size(CBX))+binCenter+xoff+SCATTER*(rand(size(CBX))-.5),subject(k).R2BX,'.','color',green,'MarkerSize',DOTSIZE)
        plot(zeros(size(CX))+binCenter+xoff+SCATTER*(rand(size(CX))-.5),subject(k).R2X,'.','color',black,'MarkerSize',DOTSIZE)
        plot(zeros(size(CY))+binCenter+xoff+SCATTER*(rand(size(CY))-.5),subject(k).R2Y,'.','color',red,'MarkerSize',DOTSIZE)

        lists(k,R).x=Xlist(:,1);
        lists(k,R).s=Xlist(:,1)*0+k;
        lists(k,R).y=Ylist(:,1);
    end
end

for k=1:NR
    %for k=1:NR
    fullX=vertcat(lists(:,k).x);
    fullY=vertcat(lists(:,k).y);
    [p,table,stats]=anova1(fullY,fullX,'off');
    kw(k).c=multcompare(stats,'display','off');

    [p,table,stats]=anova1(fullY,fullX,'off');
    mids(k,:)=stats.means;
    %     for kk=1:3
    %         F=find(fullX==kk);
    %         mids(k,kk)=median(fullY(F));
    %     end

    errors(k,:)=halfWidth(stats)';
    rangemids(k)=5*(edges(k)+edges(k+1))/2;
end
fullC=vertcat(kw.c);
calcC=[fullC(:,[1 2]) 0<=(fullC(:,4).*fullC(:,5))]

fill([rangemids rangemids(end:-1:1)],[mids(:,1)-errors(:,1); mids(end:-1:1,1)+errors(end:-1:1,1)],'w','facecolor',green,'edgecolor',green)
fill([rangemids rangemids(end:-1:1)],[mids(:,2)-errors(:,2); mids(end:-1:1,2)+errors(end:-1:1,2)],'w','facecolor',black,'edgecolor',black)
fill([rangemids rangemids(end:-1:1)],[mids(:,3)-errors(:,3); mids(end:-1:1,3)+errors(end:-1:1,3)],'w','facecolor',red,'edgecolor',red)

xlabel('Time Post Onset of Disturbing Forces, ms')
ylabel('Max Perpendicular Distance, cm')
set(gcf,'color',[1 1 1])
set(gcf,'position',[625   250   440   490])
laprint(gcf,'figures/fig4raw','scalefonts','off','asonscreen','on')



