function batchFig3Old(S)

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
yoff=[3 1; 2 0]*.13;

figure(37)
clf
hold on

load(['../Data/Data_pulse/pulse',num2str(S),'W.mat'])
load(['../Data/Data_pulse/pulse',num2str(S),'Y.mat'])

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

margin=.02;
set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
set(gcf,'color',[1 1 1])

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

%matlabfrag('figures/fig3raw')
%laprint(gcf,'figures/fig3raw','width',15,'scalefonts','off','factor',1)

print(['figures/fig3Old_',num2str(S)],'-dtiff','-r300')