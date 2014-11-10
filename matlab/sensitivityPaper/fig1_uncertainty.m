clc
clear all

qscale=.001;
SKIP=2;
yoff=.04;

gray=.5*[1 1 1];
green=[.1 .15 .7];
red=[1 .5 .1];
pink=[1 .8 .8];
blue=[.8 .8 1];
black=[0 0 0];
offset=[0 -.06 -.1];

setGlobals(paramsPopulator('burdet'))

figure(1)
clf
hold on
batches=[7 8 2];

coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);

for k=1:3
    load(['BATCH',num2str(batches(k)),'.mat'])
    
    rp=randperm(1000);
    for kk=1:50
        plot(simA(rp(kk)).y(:,1),simA(rp(kk)).y(:,2)+offset(k),'Color',pink,'linewidth',.5)
    end
    
    y=extract(t,xvaf,'reflex');
    
    tcalc=t;
    tcalc(t>=.7)=.7;
    x=minjerk(coeff,tcalc)';
    
    plot(xvaf(:,1),xvaf(:,2)+offset(k),'-','Color',black,'Linewidth',2)
    X=[xvaf(1:SKIP:end,1) xvaf(1:SKIP:end,2)+offset(k)];
    Y=X+qscale*xvaf(1:SKIP:end,7:8);
    arrow(X,Y,gray,.3)
    plot(x(1:SKIP:end,1),x(1:SKIP:end,2)+offset(k),'o-','Color',green,'markerfacecolor',green)
    plot(y(1:SKIP:end,1),y(1:SKIP:end,2)+offset(k),'.-','Color',red,'markersize',12)
    
    
end

margin=.02;
set(gcf,'position',[250 300 625 305])
set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
set(gcf,'color',[1 1 1])

axis equal
axis off

latexscale=1.2;

text(0,.54,'Early Pulse','horizontalalignment','left','Verticalalignment','top','color','k')
text(0,.48,'Late Pulse','horizontalalignment','left','Verticalalignment','top','color','k')
text(0,.4325,'Filtered Gaussian','horizontalalignment','left','Verticalalignment','top','color','k')

%annotate(h)
colors=[gray;black;green;red;red];
labs={'Force Disturbance','Hand Trajectory','Desired Hand Trajectory','Extracted Desired Hand Trajectory','Extracted with Parameter Error'};
p=[.078,.4455;.1263,.5312;.1301,.4003;.145,.40;.12,.4932];
d=[1,-1;1,1;1,-1;1,-1;1,0];
al=.012;
alength=[al;al;.8*al;1.7*al;.5*al];

h=annotate(p,d,labs,colors,alength);

plot([0 .01],.496*[1 1],'k','linewidth',3)
text(0.005,.4955,'1 cm','Horizontalalignment','center','Verticalalignment','top')
arrow([.005 .505],[.005 .505]+[0 qscale*10],gray,.3)
text(.00555,.5055,'10 N','rotation',90,'Horizontalalignment','center','Verticalalignment','bottom')

set(0,'defaulttextinterpreter','none')

matlabfrag('../figures/fig1raw');
%laprint(gcf,'../figures/fig1raw','width',15,'scalefonts','off','factor',1)

%set(0,'defaulttextinterpreter','latex')
%print ../figures/fig1.eps -depsc