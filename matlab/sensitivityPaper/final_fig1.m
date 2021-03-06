clc
clear all

if ~exist('finalfig1data.mat','file')|1
    setGlobals(paramsPopulator('burdet'))
    
    load BATCH2.mat

    t=0:.005:2;
    coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
    tcalc=t;
    tcalc(t>=.7)=.7;
    [x,v,a]=minjerk(coeff,tcalc);
    x=x';
    v=v';
    a=a';
    fnhxva=[x v a];
    wh.xvaf=xvaf;
    wh.x=x;
    wh.xsim=xvaf(:,[1 2]);
    wh.yex=extract(t,xvaf,'reflex');
    
    f=zeros(length(t),2);
    %early
    first=find(x(:,1)>=.015,1,'first');
    last=find(t>=(t(first)+.15),1,'first');
    fearly=f;
    fearly(first:last,2)=15;
    
    xsim=forwardSim(paramsPopulator,t,[fnhxva fearly]);
    yex=extract(t,[xsim fearly],'reflex');
    epul.xvaf=[xsim fearly];
    epul.x=x;
    epul.xsim=xsim;
    epul.yex=yex;
    
    f=zeros(length(t),2);
    %early
    first=find(x(:,1)>=.075,1,'first');
    last=find(t>=(t(first)+.15),1,'first');
    flate=f;
    flate(first:last,2)=15;
    
    xsim=forwardSim(paramsPopulator,t,[fnhxva flate]);
    yex=extract(t,[xsim flate],'reflex');
    lpul.xvaf=[xsim flate];
    lpul.x=x;
    lpul.xsim=xsim;
    lpul.yex=yex;
    
    save('finalfig1data.mat','wh','epul','lpul')
else
    load finalfig1data.mat
end

qscale=.001;
SKIP=2;
yoff=.04;

figure(1)
clf
hold on

gray=.5*[1 1 1];
green=[.1 .7 .3];
red=[1 .3 .3];
blue=[.3 .3 1];
black=[0 0 0];

plotme(1)=epul;
plotme(2)=lpul;
plotme(3)=wh;

offset=[0 -.05 -.07];

for k=1:2
plot(plotme(k).xsim(:,1),plotme(k).xsim(:,2)+offset(k),'-','Color',black,'Linewidth',2)
%quiver(plotme(k).xsim(1:SKIP:end,1),plotme(k).xsim(1:SKIP:end,2)+offset(k),qscale*plotme(k).xvaf(1:SKIP:end,7),qscale*plotme(k).xvaf(1:SKIP:end,8),0,'Color',gray)
X=[plotme(k).xsim(1:SKIP:end,1) plotme(k).xsim(1:SKIP:end,2)+offset(k)];
Y=X+qscale*plotme(k).xvaf(1:SKIP:end,7:8);
arrow(X,Y,gray,.3)
plot(plotme(k).x(1:SKIP:end,1),plotme(k).x(1:SKIP:end,2)+offset(k),'o-','Color',green,'markerfacecolor',green)
plot(plotme(k).yex(1:SKIP:end,1),plotme(k).yex(1:SKIP:end,2)+offset(k),'.-','Color',red,'markersize',12)
end


margin=.02;
set(gcf,'position',[250 300 625 305])
set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
set(gcf,'color',[1 1 1])

axis equal
axis off

latexscale=1.2;

x=.077;
y=.455;
l=.01*sqrt(2)/2*[-1 1];
arrow([x y]+latexscale*l,[x y],gray,.3,2)
text(x+l(1), y+l(2),'Force Disturbance','horizontalalignment','right','Verticalalignment','bottom','color',gray)

x=.1185;
y=.53;
l=-.01*sqrt(2)/2*[1 1];
arrow([x y]+latexscale*l,[x y],black,.3,2)
text(x+l(1), y+l(2),'Hand Trajectory','horizontalalignment','right','Verticalalignment','top','color',black)

x=.106;
y=.501;
l=.01*sqrt(2)/2*[-1 1];
arrow([x y]+latexscale*l,[x y],green,.3,2)
text(x+l(1), y+l(2),'Desired Hand Trajectory','horizontalalignment','right','Verticalalignment','bottom','color',green)

x=.1005;
y=.497;
l=-.01*sqrt(2)/2*[1 1];
arrow([x y]+latexscale*l,[x y],red,.3,2)
text(x+l(1), y+l(2),'Extracted Desired Hand Trajectory','horizontalalignment','right','Verticalalignment','top','color',red)


plot([0 .01],.486*[1 1],'k','linewidth',3)
text(0,.4855,'1 cm','Horizontalalignment','left','Verticalalignment','top')
arrow([.001 .505],[.001 .505]+[0 qscale*10],gray,.3)
text(.0005,.51,'10 N','rotation',90,'Horizontalalignment','center','Verticalalignment','bottom')

set(0,'defaulttextinterpreter','none')

matlabfrag('../figures/fig1raw');
%laprint(gcf,'../figures/fig1raw','width',15,'scalefonts','off','factor',1)

set(0,'defaulttextinterpreter','latex')
print ../figures/fig1.eps -depsc