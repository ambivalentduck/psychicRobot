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
edges=[0 linspace(0,200,NR+1)]; %sampling rate is 200 Hz = .005 s samples

black=[0 0 0];
red=[1 .3 .3];
green=[.1 .7 .3];
width=12;
height=1.618*width;
lmargin=.5;
bmargin=0;
figmargin=.4;

figure(1)
clf
hold on
set(gcf,'color',[1 1 1])
set(gcf,'units','centimeters')
set(gcf,'position',[4,8,figmargin+lmargin+width,figmargin+bmargin+height])

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
        cleanStruct(t).R2BX=getR2(trials(f(t)).x(inds,:));
        cleanStruct(t).R2BY=getR2(undisturbed(f(t)).y);
    end
    subject(k).R2BX=vertcat(cleanStruct.R2BX);
    subject(k).R2BY=vertcat(cleanStruct.R2BY);

    clear triStruct

    for R=1:NR+1
        RANGE=edges(R):edges(R+1);
        
        for t=1:length(trials)
            if ~sum(trials(t).disturbcat==[1:4]) %Fast-ish skipping instead of indexing
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
        subject(k).R2X=subject(k).R2X*sign(subject(k).R2X(3));
        subject(k).R2Y=vertcat(triStruct.R2Y);

        CBX=1*ones(size(subject(k).R2BX));
        CBY=4*ones(size(subject(k).R2BY));
        CX=2*ones(size(subject(k).R2X));
        CY=3*ones(size(subject(k).R2Y));

        %Xlist=[CBX;CBY;CX;CY];
        %Ylist=[subject(k).R2BX;subject(k).R2BY;subject(k).R2X;subject(k).R2Y];
        Xlist=[CBX;CX;CY];
        Ylist=[subject(k).R2BX;subject(k).R2X;subject(k).R2Y];

        DOTSIZE=3;
        LINEWIDTH=5;
        %binCenter=5*(edges(R)+edges(R+1))/2;
        binCenter=5*edges(R+1);
        drawnWidth=5*.35*(edges(end)-edges(end-1));
        SCATTER=.2*drawnWidth;
        NLines=4;
        xoff=drawnWidth/2*(2/(NLines-1)*(k-1)-1);

%         pctiles = prctile(subject(k).R2BX,[25;50;75]);
%         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',green,'LineWidth',LINEWIDTH)
%         pctiles = prctile(subject(k).R2X,[25;50;75]);
%         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',black,'LineWidth',LINEWIDTH)
%         pctiles = prctile(subject(k).R2Y,[25;50;75]);
%         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',red,'LineWidth',LINEWIDTH)

        %plot(zeros(size(CBX))+binCenter+xoff+SCATTER*(rand(size(CBX))-.5),subject(k).R2BX,'.','color',green,'MarkerSize',DOTSIZE)
        plot(zeros(size(CX))+binCenter+xoff+SCATTER*(rand(size(CX))-.5),subject(k).R2X,'.','color',black,'MarkerSize',DOTSIZE)
        plot(zeros(size(CY))+binCenter+xoff+SCATTER*(rand(size(CY))-.5),subject(k).R2Y,'.','color',red,'MarkerSize',DOTSIZE)

        lists(k,R).x=Xlist(:,1);
        lists(k,R).s=Xlist(:,1)*0+k;
        lists(k,R).y=Ylist(:,1);
    end
end

NR=200; %200
edges=[0 linspace(0,200,NR+1)]; %sampling rate is 200 Hz = .005 s samples

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
        X=trials(f(t)).x(inds,2)-.5;
        [ma,indX]=max(abs(X));
        cleanStruct(t).R2BX=X(indX)*100;
        cleanStruct(t).R2BY=getR2(undisturbed(f(t)).y);
    end
    subject(k).R2BX=vertcat(cleanStruct.R2BX);
    subject(k).R2BXP=subject(k).R2BX(subject(k).R2BX>0);
    subject(k).R2BXM=subject(k).R2BX(subject(k).R2BX<=0);
    subject(k).R2BY=vertcat(cleanStruct.R2BY);

    clear triStruct

    for R=1:NR+1
        RANGE=edges(R):edges(R+1);
        
        for t=1:length(trials)
            if ~sum(trials(t).disturbcat==[1:4]) %Fast-ish skipping instead of indexing
                continue
            end
            inds=trialInfo(t).forceinds(1)+RANGE;
            onset=find(vecmag(trials(t).v)>.05,1,'first');
            start=max(onset-35,1);
            indsy=max(1,trialInfo(t).forceinds(1)-start)+RANGE;
            
            sPeak=sign(trials(t).x(trialInfo(t).forceinds(1)+50,2)-.5);

            try
                triStruct(t).R2X=sPeak*getR2(trials(t).x(inds,:));
                %triStruct(t).R2Y=sPeak*getR2(trials(t).y(indsy,:));
                triStruct(t).R2Y=getR2(trials(t).y(indsy,:));
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

        CBX=1*ones(size(subject(k).R2BXP));
        CBY=2*ones(size(subject(k).R2BXM));
        CX=3*ones(size(subject(k).R2X));
        CY=4*ones(size(subject(k).R2Y));

        Xlist=[CBX;CBY;CX;CY];
        Ylist=[subject(k).R2BXP;subject(k).R2BXM;subject(k).R2X;subject(k).R2Y];
        %Xlist=[CBX;CX;CY];
        %Ylist=[subject(k).R2BX;subject(k).R2X;subject(k).R2Y];

        DOTSIZE=3;
        LINEWIDTH=5;
        %binCenter=5*(edges(R)+edges(R+1))/2;
        binCenter=5*edges(R+1);
        drawnWidth=5*.35*(edges(end)-edges(end-1));
        SCATTER=.2*drawnWidth;
        NLines=4;
        xoff=drawnWidth/2*(2/(NLines-1)*(k-1)-1);

%         pctiles = prctile(subject(k).R2BX,[25;50;75]);
%         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',green,'LineWidth',LINEWIDTH)
%         pctiles = prctile(subject(k).R2X,[25;50;75]);
%         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',black,'LineWidth',LINEWIDTH)
%         pctiles = prctile(subject(k).R2Y,[25;50;75]);
%         plot([0 0]+binCenter+xoff,pctiles([1 3]),'-','color',red,'LineWidth',LINEWIDTH)

        %plot(zeros(size(CBX))+binCenter+xoff+SCATTER*(rand(size(CBX))-.5),subject(k).R2BX,'.','color',green,'MarkerSize',DOTSIZE)
        %plot(zeros(size(CX))+binCenter+xoff+SCATTER*(rand(size(CX))-.5),subject(k).R2X,'.','color',black,'MarkerSize',DOTSIZE)
        %plot(zeros(size(CY))+binCenter+xoff+SCATTER*(rand(size(CY))-.5),subject(k).R2Y,'.','color',red,'MarkerSize',DOTSIZE)

        lists(k,R).x=Xlist(:,1);
        lists(k,R).s=Xlist(:,1)*0+k;
        lists(k,R).y=Ylist(:,1);
    end
end

%% From here on

for k=1:NR+1
    %for k=1:NR
    fullX=vertcat(lists(:,k).x);
    fullY=vertcat(lists(:,k).y);
    [p,table,stats]=kruskalwallis(fullY,fullX,'off');
    kw(k).c=multcompare(stats,'display','off');
    
    [p,table,stats]=anova1(fullY,fullX,'off');
    mids(k,:)=stats.means;
%     for kk=1:4
%         F=find(fullX==kk);
%         mids(k,kk)=median(fullY(F));
%     end
    
    errors34(k,:)=halfWidth(stats)';
    errors12(k,:)=halfWidth(stats,.05)';
    errors(k,:)=[errors12(k,1:2) errors34(k,3:4)];
    %rangemids(k)=5*(edges(k)+edges(k+1))/2;
    rangemids(k)=5*edges(k+1);
    kw(k).ru=rangemids(k)*[1 1 1]';
%    text(rangemids(k),-.01,num2str(rangemids(k)),'horizontalalignment','center','verticalalignment','top')
end
% fullC=vertcat(kw.c);
% calcDiff=[vertcat(kw.ru) fullC(:,[1 2]) 0<=(fullC(:,3).*fullC(:,5))]
% diffLT=[vertcat(kw.ru) fullC(:,[1 2]) fullC(:,[3 4 5])]
% diffLT(sum(diffLT(:,2:3),2)==3,:)
% diffLT(sum(diffLT(:,2:3),2)==4,:)

ALPHA=.4;
h=[fill([rangemids rangemids(end:-1:1)],[mids(:,2)-errors(:,2); mids(end:-1:1,1)+errors(end:-1:1,1)],'w','facecolor',green,'edgecolor',green);
fill([rangemids rangemids(end:-1:1)],[mids(:,3)-errors(:,3); mids(end:-1:1,3)+errors(end:-1:1,3)],'w','facecolor',black,'edgecolor',black);
fill([rangemids rangemids(end:-1:1)],[mids(:,4)-errors(:,4); mids(end:-1:1,4)+errors(end:-1:1,4)],'w','facecolor',red,'edgecolor',red)];
for k=1:3
    set(h(k),'EdgeAlpha',ALPHA,'FaceAlpha',ALPHA);
end

color=[0 0 195]/255;

plot(-40+[0,0],[0 5],'-','linewidth',3,'color',color)
text(-40,2.5,'Error, 5 cm','rotation',90,'horizontalalignment','center','verticalalignment','bottom','color',color)
set(gca,'ycolor',[1 1 1])
set(gca,'ytick',[])
set(gca,'xcolor',[1 1 1])
set(gca,'xtick',[])
set(gca,'TickLength',[0 0]);
set(gca,'linewidth',1);
xlim([-40 1020])
ylim([-13 15])
for k=0:100:1000
    text(k,-1.4,num2str(k),'horizontalalignment','center','verticalalignment','top','color',color)
end
plot([0,1000],[0 0],'-','linewidth',3,'color',color)
text(500,-10,'Time Post Onset of Disturbing Forces, ms','horizontalalignment','center','verticalalignment','middle','color',color)

set(gca,'units','centimeters')
set(gca,'position',[lmargin bmargin width height])
set(0,'defaulttextinterpreter','none')

%% Save the image
matlabfrag('figures/fig4raw','renderer','opengl','dpi',600);

%print figures/fig4.eps -depsc





