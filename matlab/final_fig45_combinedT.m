clc
clear all
close all

SUBS=1:8;
distCats=[3 4];

black=[0 0 0];
green=[.1 .15 .7];
red=[1 .5 .1];
width=12;
height=1.618*width;
lmargin=.5;
bmargin=0;
figmargin=.4;

yoffSpan=8;
plotOrder=1:8 %[1 4 2 3 5 6 7 8];

figure(1)
clf
hold on
set(gcf,'color',[1 1 1])
set(gcf,'units','centimeters')
set(gcf,'position',[4,8,figmargin+lmargin+width,figmargin+bmargin+height])

%% Compute everything we'll need later

if ~exist('fig4dotsNmeans.mat','file')|1
    for k=SUBS
        load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
        load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])
        
        %First, get plottables
        NR=10;
        edges=[0 linspace(0,200,NR+1) 200]; %sampling rate is 200 Hz = .005 s samples
        
        %Force disturbed, deviation hand and extracted
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
                    triStruct(t).dX=sPeak*getR2(trials(t).x(inds,:));
                    triStruct(t).dY=sPeak*getR2(trials(t).y(indsy,:));
                catch
                    st1=size(trials(t).x,1);
                    triStruct(t).dX=sPeak*getR2(trials(t).x(inds(1):st1,:));
                    st1=size(trials(t).y,1);
                    triStruct(t).dY=sPeak*getR2(trials(t).y(indsy(1):st1,:));
                end
                
            end
            subjectplot(k,R).dX=vertcat(triStruct.dX);
            subjectplot(k,R).dY=vertcat(triStruct.dY);
        end
        
        %Second, get t-testables
        NR=200; %200
        edges=[0 linspace(0,200,NR+1)]; %sampling rate is 200 Hz = .005 s samples
        
        %Trials not near a disturbance, deviation baseline hand
        clear cleanStruct
        f=find([trialInfo.clean]);
        for t=2:length(f) %Trial 1 isn't really a trial
            X=trials(f(t)).x(:,2)-.5;
            [ma,indX]=max(abs(X));
            cleanStruct(t).dBX=X(indX)*100;
        end
        subjectbase(k).dBX=vertcat(cleanStruct.dBX);
        
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
                    triStruct(t).dX=sPeak*getR2(trials(t).x(inds,:));
                    triStruct(t).dY=sPeak*getR2(trials(t).y(indsy,:));
                catch %Very rarely, you need to account for indexing issues
                    thishappened=1
                    st1=size(trials(t).x,1);
                    triStruct(t).dX=sPeak*getR2(trials(t).x(inds(1):st1,:));
                    st1=size(trials(t).y,1);
                    triStruct(t).dY=sPeak*getR2(trials(t).y(indsy(1):st1,:));
                end
            end
            subjecttest(k,R).dX=vertcat(triStruct.dX);
            subjecttest(k,R).dY=vertcat(triStruct.dY);
        end
        
    end
    save('fig4dotsNmeans.mat','subjectplot','subjecttest','subjectbase');
else
    load('fig4dotsNmeans.mat')
end

%% Do the plotting
NR=10;
edges=[0 linspace(0,200,NR+1) 200]; %sampling rate is 200 Hz = .005 s samples

for k=SUBS
    for R=1:NR+1
        sX=size(subjectplot(k,R).dX);
        sY=size(subjectplot(k,R).dY);
        
        DOTSIZE=3;
        LINEWIDTH=5;
        binCenter=5/2*(edges(R)+edges(R+1));
        drawnWidth=5; %2*.35*(edges(end-1)-edges(end-2));
        SCATTER=5; %*drawnWidth;
        NLines=2;
        xoff=10*mod(plotOrder(k)-1,2)-5; %drawnWidth/2*(2/(NLines-1)*(k-1)-1);
        yoff=(plotOrder(k)-1)*yoffSpan;
        
        plot(zeros(sX)+binCenter+xoff+SCATTER*(rand(sX)-.5),subjectplot(k,R).dX+yoff,'.','color',black,'MarkerSize',DOTSIZE)
        plot(zeros(sY)+binCenter+xoff+SCATTER*(rand(sY)-.5),subjectplot(k,R).dY+yoff,'.','color',red,'MarkerSize',DOTSIZE)
    end
end




%% From here on

color=[0 0 195]/255;
NR=200; %200
        edges=[0 linspace(0,200,NR+1)]; %sampling rate is 200 Hz = .005 s samples

for S=SUBS
    S
    %subplot(4,1,S)
    onset=-inf;
    handonset=-inf;
    %RP=randperm(length(subject(S).R2BX));
    %P=subject(S).R2BX(RP(1:most));
    %P=sort(subjectbase(S).dBX,1,'descend');
    P=abs(subjectbase(S).dBX);
    %P=P(70+(1:most));
    %RP=randperm(length(P));
    %P=P(RP(1:20));
    %M=subject(S).R2BXM(RP(1:most));
    mids=zeros(NR+1,3);
    lowers=zeros(NR+1,3);
    uppers=zeros(NR+1,3);
    mP=mean(P);
    uppers(:,1)=ones(NR+1,1)*mP;
    lowers(:,1)=-ones(NR+1,1)*mP;
    P=P-mP; %Simplifies plotting confidence intervals
    yoff=(plotOrder(S)-1)*yoffSpan;
    
    for k=1:NR+1
        x=subjecttest(S,k).dX;
        x=x(~isnan(x));
        mx=mean(x);
        sem95=1.684*std(x)/sqrt(length(x));
        mids(k,2)=mx;
        uppers(k,2)=mx+sem95;
        lowers(k,2)=mx-sem95;
        
        if mx>0
            [xn(S,k),pr,cix]=ttest2(P,x-mP,.05,'left','unequal');
            %uppers(k,2)=mx;
            %lowers(k,2)=max(mx+ci(2),0);
        else
            [xn(S,k),pr,cix]=ttest2(P,x+mP,.05,'right','unequal');
            %uppers(k,2)=mx;
            %lowers(k,2)=min(0,mx+sem95);
        end
            
        y=subjecttest(S,k).dY;
        my=mean(y);
        sem95=1.684*std(y)/sqrt(length(y));
        mids(k,3)=my;
        uppers(k,3)=my+sem95;
        lowers(k,3)=my-sem95;
        if my>0
            [yn(S,k),pr,ciy]=ttest2(P,y-mP,.05,'left','unequal');
            %uppers(k,3)=my;
            %lowers(k,3)=min(0,my-sem95);
        else
            [yn(S,k),pr,ciy]=ttest2(P,y+mP,.05,'right','unequal');
            %uppers(k,3)=my;
            %lowers(k,3)=max(my+sem95,0);
        end
        
        lscale=1.5;
        rangemids(k)=5*edges(k+1);
        if onset<0
            if yn(S,k)
                onset=rangemids(k)
                plot(rangemids(k)+[0 0],lscale*mP*[-1 1]+yoff,'color',[1 0 0],'linewidth',1.5)
                %plot(rangemids(k),.5+lscale*mP+yoff,'o','markerfacecolor',[1 0 0],'markersize',10)
            end
        end
        if handonset<0
            if xn(S,k)
                handonset=rangemids(k)
                plot(rangemids(k)+[0 0],lscale*mP*[-1 1]+yoff,'color',[0 1 0],'linewidth',1.5)
                %plot(rangemids(k),.5+lscale*mP+yoff,'o','markerfacecolor',[0 1 0],'markersize',10)
            end
        end
    end
    differences(S)=onset-handonset

    plot([0,1000],yoff+[0 0],'-','linewidth',1,'color',color)
    
    rangemids(1)=-5;
    rangemids(end)=1005;
    
    ALPHA=.4;
    h=[fill([rangemids rangemids(end:-1:1)],[lowers(:,1); uppers(end:-1:1,1)]+yoff,'w','facecolor',green,'edgecolor',green);
        fill([rangemids rangemids(end:-1:1)],[lowers(:,2); uppers(end:-1:1,2)]+yoff,'w','facecolor',black,'edgecolor',black);
        fill([rangemids rangemids(end:-1:1)],[lowers(:,3); uppers(end:-1:1,3)]+yoff,'w','facecolor',red,'edgecolor',red)];
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
ylim([-10 60])
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





