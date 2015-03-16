clc
clear all

figure(37)
clf
subplot(5,1,1:3)
hold on
specialestsub=5;

%% Fig3 stuff

lw=1.5;

SKIP=9;
qscale=.002;
fade=.7;

gray=.5*[1 1 1];
green=[.1 .15 .7];
red=[1 .5 .1];
pink=[1 .8 .8];
darkpink=.8*[.8 1 .8];
blue=[.8 .8 1];
black=[0 0 0];

xoff=[.15 .15;0 0];
yoff=[3 1; 2 0]*.1;

S=specialestsub;
load(['../Data/Data_pulse/pulse',num2str(S),'W.mat'])
load(['../Data/Data_pulse/pulse',num2str(S),'Y.mat'])
%load(['../Data/Data_pulse/pulse',num2str(S),'U.mat'])

starts=[trialInfo.startcat];
ends=[trialInfo.endcat];
dcat=[trials.disturbcat];
clean=[trialInfo.clean];

f=find(clean);
for c=2:length(f)
    if rand<.8
        continue
    end
    k=f(c);
    X=[trials(k).x(1:SKIP:end,1)-means(starts(k)),trials(k).x(1:SKIP:end,2)-.5];
    rlength=abs(starts(k)-ends(k));
    rdir=sign(ends(k)-starts(k));
    if rdir<0
        X=-X;
    end
    for earlylate=1:2
        xo=xoff(rlength,earlylate);
        yo=yoff(rlength,earlylate);
        plot(X(:,1)+xo,X(:,2)+yo,'color',green,'linewidth',.5)
    end
end

for k=2:length(trials)
    %First decide length class
    rlength=abs(starts(k)-ends(k));
    rdir=sign(ends(k)-starts(k));
    switch dcat(k)
        case {1 2}
            earlylate=1;
        case {3 4}
            earlylate=2;
        otherwise
            %Not a pulse: white or undisturbed
            continue
    end
    X=[trials(k).x(:,1)-means(starts(k)),trials(k).x(:,2)-.5];
    Y=[trials(k).y(:,1)-means(starts(k)),trials(k).y(:,2)-.5];
    F=trials(k).f(1:SKIP:end,:);
    if rdir<0
        continue
        X=-X;
        Y=-Y;
        F=-F;
    else
        %continue
    end
    if mod(dcat(k),2)
        X(:,2)=-X(:,2);
        Y(:,2)=-Y(:,2);
        F(:,2)=-F(:,2);
    end

    xo=xoff(rlength,earlylate);
    yo=yoff(rlength,earlylate);
    plot(X(:,1)+xo,X(:,2)+yo,'color',black,'linewidth',lw)
    plot(Y(:,1)+xo,Y(:,2)+yo,'color',red,'linewidth',lw)
    XF=[X(1:SKIP:end,1)+xo X(1:SKIP:end,2)+yo];
    arrow(XF,XF+qscale*F,gray,.3,lw);
end

%margin=.02;
%set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
%set(gcf,'color',[1 1 1])

axis equal

%annotate(h);
p=[.18 .055;.2092 -.01458;.28 .3619;.2 .302];
d=[1 0; 1 1;1 0;-1 0];
alength=[.03,.02,.03,.03];
colors=[gray; green; red; black];
labs{1}='Force Disturbance';
labs{2}='Undisturbed Hand Trajectory';
labs{3}='Extracted Desired Hand Trajectory';
labs{4}='Hand Trajectory';

h=annotate(p,d,labs,colors,alength);

lshift=.01;
plot([.01 .04]-lshift,-.01*[1 1],'k','linewidth',3)
text(.025-lshift,-.011,'3 cm','horizontalalignment','center','verticalalignment','top')

A=.008*[1 1];
A(1)=A(1)-lshift;
arrow(A,A+[0 qscale*10],gray,.3,2)
text(.006-lshift,(.008+qscale*5),'10 N','rotation',90,'Horizontalalignment','center','Verticalalignment','bottom')

f=findobj('Type','text');
set(f,'fontsize',10)


text(0,.46,'Early Pulse Disturbance','horizontalalignment','left','verticalalignment','top','fontsize',12,'fontweight','bold');
%text(0,.46,'Early Pulse Disturbance','horizontalalignment','left','verticalalignment','top');
text(0,.2,'Late Pulse Disturbance','horizontalalignment','left','verticalalignment','top','fontsize',12,'fontweight','bold');


axis off

set(gcf,'units','centimeters')
%set(gcf,'position',[3 3 18 18])

%% Fig4 stuff

load('fig4dotsNmeans.mat')
subplot(5,1,4)
hold on
NR=10;
DOTSIZE=3;
LINEWIDTH=5;
SCATTER=5; %*drawnWidth;
NLines=2;

edges=[0 linspace(0,200,NR+1) 200]; %sampling rate is 200 Hz = .005 s samples


for k=specialestsub
    for R=1:NR+1
        sX=size(subjectplot(k,R).dX);
        sY=size(subjectplot(k,R).dY);

        binCenter=5/2*(edges(R)+edges(R+1));
        drawnWidth=5; %2*.35*(edges(end-1)-edges(end-2));

        xoff=0; %10*mod(plotOrder(k)-1,2)-5; %drawnWidth/2*(2/(NLines-1)*(k-1)-1);
        yoff=0; %(plotOrder(k)-1)*yoffSpan;

        plot(zeros(sX)+binCenter+xoff+SCATTER*(rand(sX)-.5),subjectplot(k,R).dX+yoff,'.','color',black,'MarkerSize',DOTSIZE)
        plot(zeros(sY)+binCenter+xoff+SCATTER*(rand(sY)-.5),subjectplot(k,R).dY+yoff,'.','color',red,'MarkerSize',DOTSIZE)
    end
end

color=[0 0 195]/255;
NR=200; %200
edges=[0 linspace(0,200,NR+1)]; %sampling rate is 200 Hz = .005 s samples

for S=1:8
    S
    onset=-inf;
    handonset=-inf;
    P=abs(subjectbase(S).dBX);
    mids=zeros(NR+1,3);
    lowers=zeros(NR+1,3);
    uppers=zeros(NR+1,3);
    mP=mean(P);
    uppers(:,1)=ones(NR+1,1)*mP;
    lowers(:,1)=-ones(NR+1,1)*mP;
    P=P-mP;

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
        else
            [xn(S,k),pr,cix]=ttest2(P,x+mP,.05,'right','unequal');
        end

        y=subjecttest(S,k).dY;
        my=mean(y);
        sem95=1.684*std(y)/sqrt(length(y));
        mids(k,3)=my;
        uppers(k,3)=my+sem95;
        lowers(k,3)=my-sem95;
        if my>0
            [yn(S,k),pr,ciy]=ttest2(P,y-mP,.05,'left','unequal');
        else
            [yn(S,k),pr,ciy]=ttest2(P,y+mP,.05,'right','unequal');
        end

        lscale=1.5;
        rangemids(k)=5*edges(k+1);

        if onset<0
            if yn(S,k)
                onset=rangemids(k)
            end
        end
        if handonset<0
            if xn(S,k)
                handonset=rangemids(k)
            end
        end
        if S==specialestsub
            rangemids(1)=-5;
            rangemids(end)=1005;
            plot(onset+[0 0],lscale*mP*[-1 1]+yoff,'color',[1 0 0],'linewidth',1.5)
            plot(rangemids(k)+[0 0],lscale*mP*[-1 1]+yoff,'color',[0 1 0],'linewidth',1.5)
            ALPHA=.4;
            h=[fill([rangemids rangemids(end:-1:1)],[lowers(:,1); uppers(end:-1:1,1)]+yoff,'w','facecolor',green,'edgecolor',green);
                fill([rangemids rangemids(end:-1:1)],[lowers(:,2); uppers(end:-1:1,2)]+yoff,'w','facecolor',black,'edgecolor',black);
                fill([rangemids rangemids(end:-1:1)],[lowers(:,3); uppers(end:-1:1,3)]+yoff,'w','facecolor',red,'edgecolor',red)];
            for k=1:3
                set(h(k),'EdgeAlpha',ALPHA,'FaceAlpha',ALPHA);
            end
        end
    end
    onsets(S)=onset;
    handonsets(S)=handonset;
    differences(S)=onset-handonset

    %plot([0,1000],yoff+[0 0],'-','linewidth',1,'color',color)
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
%ylim([-10 60])
for k=0:100:1000
    text(k,-1.4,num2str(k),'horizontalalignment','center','verticalalignment','top','color',color)
end

text(500,-6,'Time Post Onset of Disturbing Forces, ms','horizontalalignment','center','verticalalignment','middle','color',color)

subplot(5,1,5)
hold on
for k=1:8
    plot(onsets(k),k,'.','color',[1 0 0])
    plot(handonsets(k),k,'.','color',[0 1 0])
    if handonsets(k)<1000
        plot([onsets(k) handonsets(k)],k+[0 0],'color',black,'linewidth',1.5)
    else
        arrow
    end
end
%matlabfrag('figures/fig3raw')