clc
clear all

anecdote(1).sub=1;
anecdote(1).SE=[1 3];
anecdote(1).dcat=2;
anecdote(1).rot=0;
anecdote(1).yoff=0;
% anecdote(1).sub=1;
% anecdote(1).SE=[1 3];
% anecdote(1).dcat=2;
% anecdote(1).rot=0;
% anecdote(1).yoff=0;
anecdote(2).sub=1;
anecdote(2).SE=[3 1];
anecdote(2).dcat=4;
anecdote(2).rot=pi;
anecdote(2).yoff=-.15;

anecdote(3).sub=1;
anecdote(3).SE=[2 3];
anecdote(3).dcat=2;
anecdote(3).rot=pi;
anecdote(3).yoff=-.3;

lw=1.5;

SKIP=4;
qscale=.002;
fade=.7;

gray=.5*[1 1 1];
green=[.1 .15 .7];
red=[1 .5 .1];
pink=[1 .8 .8];
darkpink=.8*[.8 1 .8];
blue=[.8 .8 1];
black=[0 0 0];

figure(37)
clf
hold on

for A=1:length(anecdote)

    S=anecdote(A).sub;
    load(['../Data/Data_pulse/pulse',num2str(S),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(S),'Y.mat'])
    load(['../Data/Data_pulse/pulse',num2str(S),'U.mat'])

    starts=[trialInfo.startcat];
    ends=[trialInfo.endcat];
    dcat=[trials.disturbcat];
    clean=[trialInfo.clean];

    f=find((starts==anecdote(A).SE(1))&(ends==anecdote(A).SE(2))&clean);
    m=means(anecdote(A).SE(1));
    cosa=cos(anecdote(A).rot);
    sina=sin(anecdote(A).rot);
    rotmat=[cosa sina;-sina cosa];
    for c=1:length(f)
        k=f(c);
        X=[trials(k).x(:,1)-m,trials(k).x(:,2)-.5]*rotmat;
        plot(X(:,1),X(:,2)+anecdote(A).yoff,'color',green,'linewidth',.1)
        %plot(undisturbed(k).y(:,1),undisturbed(k).y(:,2),'-','linewidth',.0000001,'color',[1 fade fade])
    end

    %Look at each reach cat.
    f=find((starts==anecdote(A).SE(1))&(ends==anecdote(A).SE(2))&(dcat==anecdote(A).dcat));
    for c=1:length(f)
        k=f(c);
        x=[trials(k).x(:,1)-m,trials(k).x(:,2)-.5]*rotmat;
        plot(x(:,1),x(:,2)+anecdote(A).yoff,'k','linewidth',lw)
        X=[x(1:SKIP:end,1) x(1:SKIP:end,2)+anecdote(A).yoff];
        F=trials(k).f(1:SKIP:end,:)*rotmat;
        arrow(X,X+qscale*F,gray,.3,lw);
        Y=[trials(k).y(:,1)-m trials(k).y(:,2)-.5]*rotmat;
        plot(Y(:,1),Y(:,2)+anecdote(A).yoff,'color',red,'linewidth',lw)
    end
end
axis equal

margin=.02;
set(gca,'position',[margin margin 1-margin 1-margin],'units','normalized')
set(gcf,'color',[1 1 1])

axis equal

%annotate(h);
p=[.17 -.105;.17 .006;.15 -.152;.123 .025];
d=[1 -1; -1 -1;1 1;-1 -1];
colors=[gray; green; red; 0 0 0];
labs{1}='Force Disturbance';
labs{2}='Undisturbed Hand Trajectory';
labs{3}='Extracted Desired Hand Trajectory';
labs{4}='Hand Trajectory';

h=annotate(p,d,labs,colors,.015);


plot([.01 .04],-.01*[1 1],'k','linewidth',3)
text(.025,-.011,'3 cm','horizontalalignment','center','verticalalignment','top')

A=.008*[1 1];
arrow(A,A+[0 qscale*10],gray,.3,2)
text(.006,(.008+qscale*5),'10 N','rotation',90,'Horizontalalignment','center','Verticalalignment','bottom')

axis off

set(gcf,'units','centimeters')
get(gcf,'position')

matlabfrag('figures/fig3raw')
%laprint(gcf,'figures/fig3raw','width',15,'scalefonts','off','factor',1)

print figures/fig3 -depsc