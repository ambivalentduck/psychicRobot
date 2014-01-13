clc
clear all

anecd=[6 7 3];

if ~exist('finalfig1data.mat','file')||1
    load BATCH3.mat
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
    
    
    f=zeros(length(t),2);
    f((t>=.1)&(t<=.15),2)=15;
    
    xsim,xsimSM]=forwardSim(paramsPopulator,t,xvaf);
    yex=extract(t,[xsim f],'reflex');
    pul.xvaf=xvaf;
    pul.x=x;
    pul.xsim=xsim;
    pul.yex=yex;
    

    
    save('finalfig1data.mat','wh','pul')
else
    load finalfig1data.mat
end
6
qscale=.001;
SKIP=2;
yoff=.04;

figure(1)
clf
hold on

plot(wh.xsim(:,1),wh.xsim(:,2),'k-')
plot(pul.xsim(:,1),pul.xsim(:,2)+yoff,'k-')
quiver(wh.xsim(1:SKIP:end,1),wh.xsim(1:SKIP:end,2),qscale*wh.xvaf(1:SKIP:end,7),qscale*wh.xvaf(1:SKIP:end,8),0,'Color',.6*[1 1 1])
quiver(pul.xsim(1:SKIP:end,1),pul.xsim(1:SKIP:end,2)+yoff,qscale*pul.xvaf(1:SKIP:end,7),qscale*pul.xvaf(1:SKIP:end,8),0,'Color',.6*[1 1 1])

exgray=.4;
ingray=.05;
plot(wh.x(:,1),wh.x(:,2),'-','color',ingray*[1 1 1],'linewidth',6)
plot(pul.x(:,1),pul.x(:,2)+yoff,'-','color',ingray*[1 1 1],'linewidth',6)
plot(wh.yex(:,1),wh.yex(:,2),'.','color',exgray*[1 1 1],'linewidth',1)
plot(pul.yex(:,1),pul.yex(:,2)+yoff,'.','color',exgray*[1 1 1],'linewidth',1)

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


plot([0 .15],.486*[1 1],'k','linewidth',3)
text(0,.4855,'15 cm','Horizontalalignment','left','Verticalalignment','top')
quiver([.0005],[.505],[0],qscale*[10],0,'Color',.6*[1 1 1],'linewidth',1)
text(.0005,.51,'10 N','rotation',90,'Horizontalalignment','center','Verticalalignment','bottom')

set(0,'defaulttextinterpreter','none')

laprint(gcf,'fig1raw','width',15,'scalefonts','off','factor',1)