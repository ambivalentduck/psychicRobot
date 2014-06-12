clc
clear all
close all

%Set ranges
%For each range, show the bar plots with figure window titled
%    set(gcf,'name',range)
%Grab the mean and confidence interval from c=multcompare with output off
%Plot mean and confidence interval against range
%call that finalfig5

SUBS=1:4;
distCats=[3 4];

NR=10;
edges=[0 linspace(0,200,NR+1) 200]; %sampling rate is 200 Hz = .005 s samples

black=[0 0 0];
red=[1 .3 .3];
green=[.1 .7 .3];
width=12;
height=1.618*width;
lmargin=.5;
bmargin=0;
figmargin=.4;

yoffSpan=8;
plotOrder=[1 4 2 3];

figure(1)
clf
hold on
set(gcf,'color',[1 1 1])
set(gcf,'units','centimeters')
set(gcf,'position',[4,8,figmargin+lmargin+width,figmargin+bmargin+height])

for k=SUBS
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])

    clear cleanStruct
    f=find([trialInfo.clean]);
    for t=2:length(f)
        X=trials(f(t)).x(:,2)-.5;
        [ma,indX]=max(abs(X));
        cleanStruct(t).R2BX=X(indX)*100;
    end
    subject(k).R2BX=vertcat(cleanStruct.R2BX);

    clear triStruct

    for R=1:NR+2
        RANGE=edges(R):edges(R+1);

        for t=1:length(trials)
            if ~sum(trials(t).disturbcat==[distCats]) %Fast-ish skipping instead of indexing
                continue
            end
            inds=trialInfo(t).forceinds(1)+RANGE;
            onset=find(vecmag(trials(t).v)>.05,1,'first');
            start=max(onset-35,1);
            indsy=max(1,trialInfo(t).forceinds(1)-start)+RANGE;

            sPeak=sign(trials(t).x(trialInfo(t).forceinds(1)+50,2)-.5);

            try
                triStruct(t).R2X=sPeak*getR2(trials(t).x(inds,:));
                triStruct(t).R2Y=sPeak*getR2(trials(t).y(indsy,:));
            catch
                st1=size(trials(t).x,1);
                triStruct(t).R2X=sPeak*getR2(trials(t).x(inds(1):st1,:));
                st1=size(trials(t).x,1);
                triStruct(t).R2Y=sPeak*getR2(trials(t).y(indsy(1):st1,:));
            end

        end
        subject(k).R2X=vertcat(triStruct.R2X);
        subject(k).R2X=subject(k).R2X;
        subject(k).R2Y=vertcat(triStruct.R2Y);

        CBX=1*ones(size(subject(k).R2BX));
        CX=2*ones(size(subject(k).R2X));
        CY=3*ones(size(subject(k).R2Y));

        DOTSIZE=3;
        LINEWIDTH=5;
        binCenter=5/2*(edges(R)+edges(R+1));
        drawnWidth=5; %2*.35*(edges(end-1)-edges(end-2));
        SCATTER=5; %*drawnWidth;
        NLines=2;
        xoff=10*mod(plotOrder(k)-1,2)-5; %drawnWidth/2*(2/(NLines-1)*(k-1)-1);
        yoff=(plotOrder(k)-1)*yoffSpan;

        %plot(zeros(size(CBX))+binCenter+xoff+SCATTER*(rand(size(CBX))-.5),subject(k).R2BX+yoff,'.','color',green,'MarkerSize',DOTSIZE)
        plot(zeros(size(CX))+binCenter+xoff+SCATTER*(rand(size(CX))-.5),subject(k).R2X+yoff,'.','color',black,'MarkerSize',DOTSIZE)
        plot(zeros(size(CY))+binCenter+xoff+SCATTER*(rand(size(CY))-.5),subject(k).R2Y+yoff,'.','color',red,'MarkerSize',DOTSIZE)
    end
end

NR=200; %200
edges=[0 linspace(0,200,NR+1)]; %sampling rate is 200 Hz = .005 s samples

for k=SUBS
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])

    flag=0;

    clear cleanStruct
    f=find([trialInfo.clean]);
    for t=2:length(f)
        X=trials(f(t)).x(:,2)-.5;
        [ma,indX]=max(abs(X));
        cleanStruct(t).R2BX=X(indX)*100;
    end
    subject(k).R2BX=abs(vertcat(cleanStruct.R2BX));
    %subject(k).R2BXP=subject(k).R2BX(subject(k).R2BX>0);
    %subject(k).R2BXM=subject(k).R2BX(subject(k).R2BX<=0);

    clear triStruct

    for R=1:NR+1
        RANGE=edges(R):edges(R+1);

        for t=1:length(trials)
            if ~sum(trials(t).disturbcat==([distCats])) %Fast-ish skipping instead of indexing
                continue
            end
            inds=trialInfo(t).forceinds(1)+RANGE;
            onset=find(vecmag(trials(t).v)>.05,1,'first');
            start=max(onset-35,1);
            indsy=max(1,trialInfo(t).forceinds(1)-start)+RANGE;

            sPeak=sign(trials(t).x(trialInfo(t).forceinds(1)+50,2)-.5);

            try
                triStruct(t).R2X=sPeak*getR2(trials(t).x(inds,:));
                triStruct(t).R2Y=sPeak*getR2(trials(t).y(indsy,:));
            catch
                thishappened=1
                st1=size(trials(t).x,1);
                triStruct(t).R2X=sPeak*getR2(trials(t).x(inds(1):st1,:));
                st1=size(trials(t).x,1);
                triStruct(t).R2Y=sPeak*getR2(trials(t).y(indsy(1):st1,:));
            end
        end
        subject(k).R2X=vertcat(triStruct.R2X);
        subject(k).R2Y=vertcat(triStruct.R2Y);

        CX=3*ones(size(subject(k).R2X));
        CY=4*ones(size(subject(k).R2Y));

        Xlist=[CX;CY];
        Ylist=[subject(k).R2X;subject(k).R2Y];

        lists(k,R).x=Xlist(:,1);
        lists(k,R).y=Ylist(:,1);
    end
end

%% From here on

color=[0 0 195]/255;

for S=SUBS
    %subplot(4,1,S)
    onset=-inf;
    handonset=-inf;
    most=80; %max(sum(lists(S,k).x==3),sum(lists(S,k).x==4));
    %RP=randperm(length(subject(S).R2BX));
    %P=subject(S).R2BX(RP(1:most));
    P=sort(subject(S).R2BX,1,'descend');
    P=P(70+(1:most));
    %RP=randperm(length(subject(S).R2BX));
    %M=subject(S).R2BXM(RP(1:most));

    for k=1:NR+1
        fullX=[ones(most,1);vertcat(lists(S,k).x)];
        fullY=[P;vertcat(lists(S,k).y)];
        %[p,table,stats]=kruskalwallis(fullY(SeS),fullX(SeS),'off');
        %kw(k).c=multcompare(stats,'display','off');

        [p,table,stats]=anova1(fullY,fullX,'off');
        mids(k,:)=stats.means;
        stats2=stats;
        stats2.means(1)=-stats2.means(1);

        errors(k,:)=halfWidth(stats,.05)';
        rangemids(k)=5*edges(k+1);
        if onset<0
            cp=multcompare(stats,'display','off','alpha',.05); %Above one or below the other, one-sided
            cm=multcompare(stats2,'display','off','alpha',.05); %Above one or below the other, one-sided
            %Is group 4 > group 1 (row 3 negative in column 5)
            %OR is group 4 < group 2 (row 5 positive in column 3)
            if (cp(2,5)<0)||(cm(2,3)>0)
                S
                onset=rangemids(k)
            end
        end
        if handonset<0
            cp=multcompare(stats,'display','off','alpha',.05); %Above one or below the other, one-sided
            cm=multcompare(stats2,'display','off','alpha',.05); %Above one or below the other, one-sided
            %Is group 4 > group 1 (row 3 negative in column 5)
            %OR is group 4 < group 2 (row 5 positive in column 3)
            if (cp(1,5)<0)||(cm(1,3)>0)
                S
                handonset=rangemids(k)
            end
        end
    end

    yoff=(plotOrder(S)-1)*yoffSpan;
    plot([0,1000],yoff+[0 0],'-','linewidth',1,'color',color)

    rangemids(1)=-5;
    rangemids(end)=1005;

    ALPHA=.4;
    h=[fill([rangemids rangemids(end:-1:1)],[-mids(:,1)-errors(:,1); mids(end:-1:1,1)+errors(end:-1:1,1)]+yoff,'w','facecolor',green,'edgecolor',green);
        fill([rangemids rangemids(end:-1:1)],[mids(:,2)-errors(:,2); mids(end:-1:1,2)+errors(end:-1:1,2)]+yoff,'w','facecolor',black,'edgecolor',black);
        fill([rangemids rangemids(end:-1:1)],[mids(:,3)-errors(:,3); mids(end:-1:1,3)+errors(end:-1:1,3)]+yoff,'w','facecolor',red,'edgecolor',red)];
    for k=1:3
        set(h(k),'EdgeAlpha',ALPHA,'FaceAlpha',ALPHA);
    end
end

plot(-30+[0,0],[0 5],'-','linewidth',1,'color',color)
text(-30,2.5,'Error, 5 cm','rotation',90,'horizontalalignment','center','verticalalignment','bottom','color',color)
set(gca,'ycolor',[1 1 1])
set(gca,'ytick',[])
set(gca,'xcolor',[1 1 1])
set(gca,'xtick',[])
set(gca,'TickLength',[0 0]);
%set(gca,'linewidth',1);
xlim([-40 1020])
ylim([-10 40])
for k=0:100:1000
    text(k,-1.4,num2str(k),'horizontalalignment','center','verticalalignment','top','color',color)
end

text(500,-6,'Time Post Onset of Disturbing Forces, ms','horizontalalignment','center','verticalalignment','middle','color',color)

set(gca,'units','centimeters')
set(gca,'position',[lmargin bmargin width height])
set(0,'defaulttextinterpreter','none')

%% Save the image
matlabfrag('figures/fig4raw','renderer','opengl','dpi',600);

%print figures/fig4.eps -depsc





