clc
clear all

S=1;

load(['../Data/Data_pulse/pulse',num2str(S),'W.mat'])
load(['../Data/Data_pulse/pulse',num2str(S),'Y.mat'])
load(['../Data/Data_pulse/pulse',num2str(S),'U.mat'])

starts=[trialInfo.startcat];
ends=[trialInfo.endcat];
dcat=[trials.disturbcat];
clean=[trialInfo.clean];

figure(37)
clf
hold on

f=find((starts==1)&(ends==3)&clean);
for c=1:length(f)
    k=f(c);
    fade=.9;
    plot(trials(k).x(:,1),trials(k).x(:,2),'color',fade*[1 1 1],'linewidth',.0000001)
    plot(undisturbed(k).y(:,1),undisturbed(k).y(:,2),'-','linewidth',.0000001,'color',[1 fade fade])
end


%Look at each reach cat.
f=find((starts==1)&(ends==3)&(dcat==2)); %S1, D2 is best early 1->3

gray=.5*[1 1 1];
green=[.1 .7 .3];
red=[1 .3 .3];

for c=1:length(f)
    k=f(c);
    SKIP=4;
    qscale=.002;
    plot(trials(k).x(:,1),trials(k).x(:,2),'k')
    X=trials(k).x(1:SKIP:end,[1 2]);
    arrow(X,X+qscale*trials(k).f(1:SKIP:end,[1 2]),gray,.3)
    plot(trials(k).y(:,1),trials(k).y(:,2),'color',red)
end
axis equal

clc
clear all

S=4;

load(['../Data/Data_pulse/pulse',num2str(S),'W.mat'])
load(['../Data/Data_pulse/pulse',num2str(S),'Y.mat'])
load(['../Data/Data_pulse/pulse',num2str(S),'U.mat'])

starts=[trialInfo.startcat];
ends=[trialInfo.endcat];
dcat=[trials.disturbcat];
clean=[trialInfo.clean];

yoff=.14;

f=find((starts==3)&(ends==1)&clean);
for c=1:length(f)
    k=f(c);
    fade=.9;
    plot(trials(k).x(:,1),trials(k).x(:,2)+yoff,'color',fade*[1 1 1],'linewidth',.0000001)
    %plot(undisturbed(k).y(:,1),undisturbed(k).y(:,2)+yoff,'-','linewidth',.0000001,'color',[1 fade fade])
end

%Look at each reach cat.
f=find((starts==3)&(ends==1)&(dcat==4)); %S4, D4 is best late
gray=.5*[1 1 1];
green=[.1 .7 .3];
red=[1 .3 .3];

for c=1:length(f)
    k=f(c);
    SKIP=4;
    qscale=.002;
    plot(trials(k).x(:,1),trials(k).x(:,2)+yoff,'k')
    X=[trials(k).x(1:SKIP:end,1) trials(k).x(1:SKIP:end,2)+yoff];
    arrow(X,X+qscale*trials(k).f(1:SKIP:end,[1 2]),gray,.3)
    plot(trials(k).y(:,1),trials(k).y(:,2)+yoff,'color',red)
end

margin=.02;
set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
set(gcf,'color',[1 1 1])

axis equal

latexscale=1.2;

x=-.08;
y=.57;
l=.01*sqrt(2)/2*[1 -1];
arrow([x y]+latexscale*l,[x y],gray,.3,2)
text(x+l(1), y+l(2),'Force Disturbance','horizontalalignment','left','Verticalalignment','top','color',gray)


plot([-.17 -.14],.45*[1 1],'k','linewidth',3)
text(-.155,.45,'3 cm','horizontalalignment','center','verticalalignment','top')

A=[-.17,.52];
arrow(A,A+[0 qscale*10],gray,.3)

axis off
