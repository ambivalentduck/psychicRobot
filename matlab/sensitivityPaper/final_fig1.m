clc
clear all

if ~exist('finalfig1data.mat','file')
    setGlobals(paramsPopulator('burdet'))
    
    load BATCH10.mat
    wh.xvaf=xvaf;
    wh.xsim=xvaf(:,[1 2]);
    wh.yex=extract(t,xvaf,'reflex');

    t=0:.005:2;
    coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
    tcalc=t;
    tcalc(t>=.7)=.7;
    [x,v,a]=minjerk(coeff,tcalc);
    x=x';
    v=v';
    a=a';
    fnhxva=[x v a];
    wh.x=x;
    
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
green=[.3 1 .3];
red=[1 .3 .3];
blue=[.3 .3 1];

plot(wh.xsim(:,1),wh.xsim(:,2),'k-')
plot(epul.xsim(:,1),epul.xsim(:,2)+yoff,'k-')
quiver(wh.xsim(1:SKIP:end,1),wh.xsim(1:SKIP:end,2),qscale*wh.xvaf(1:SKIP:end,7),qscale*wh.xvaf(1:SKIP:end,8),0,'Color',.6*[1 1 1])
quiver(epul.xsim(1:SKIP:end,1),epul.xsim(1:SKIP:end,2)+yoff,qscale*epul.xvaf(1:SKIP:end,7),qscale*epul.xvaf(1:SKIP:end,8),0,'Color',.6*[1 1 1])


exgray=.4;
ingray=.05;
plot(wh.x(:,1),wh.x(:,2),'-','color',ingray*[1 1 1],'linewidth',6)
plot(epul.x(:,1),epul.x(:,2)+yoff,'-','color',ingray*[1 1 1],'linewidth',6)
plot(wh.yex(:,1),wh.yex(:,2),'.','color',exgray*[1 1 1],'linewidth',1)
plot(epul.yex(:,1),epul.yex(:,2)+yoff,'.','color',exgray*[1 1 1],'linewidth',1)

margin=.02;
set(gcf,'position',[250 300 625 305])
set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
set(gcf,'color',[1 1 1])

axis equal
axis off

fa=annotation('textarrow');
set(fa,'position',[0.0960    0.9279   -0.0320   -0.0656])
set(fa,'string','Force Disturbance')

ha=annotation('textarrow');
set(ha,'position',[ 0.2704    0.8787   -0.0258   -0.0670])
set(ha,'string','Simulated Hand Trajectory')

ia=annotation('textarrow');
set(ia,'position',[ 0.3824    0.7837   -0.0368   -0.0984])
set(ia,'string','Intended Hand Trajectory')

ea=annotation('textarrow');
set(ea,'position',[0.6288    0.7443   -0.0240   -0.0722])
set(ea,'string','Extracted Intended Hand Trajectory')


plot([0 .01],.486*[1 1],'k','linewidth',3)
text(0,.4855,'1 cm','Horizontalalignment','left','Verticalalignment','top')
quiver([.0005],[.505],[0],qscale*[10],0,'Color',.6*[1 1 1],'linewidth',1)
text(.0005,.51,'10 N','rotation',90,'Horizontalalignment','center','Verticalalignment','bottom')

set(0,'defaulttextinterpreter','none')

laprint(gcf,'fig1raw','width',15,'scalefonts','off','factor',1)