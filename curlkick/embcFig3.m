clc
clear all

%% Set constants and initialize

%Measurements in cm. Movements are 15 cm long

[hand,intent]=getColorScheme;
gray=.5*[1 1 1];

H=20;
W=18; %Do NOT adjust

yw=15*cosd(30);

leftMargin=0;
legendThickness=2;
legendSpacing=7;
legendp=[-yw 61];
center=[0 48];
spacing=50;
font='Arial';
fontsize=12;
interpreter='tex';
handfudge=0;
intentfudge=0;


f=figure(1);
clf
set(f,'units','inches');
set(f,'Position',[2,2,7,2.25])

hold on

tlist=(1:480)';
b=floor((tlist-1)/96)+1;

%% Draw Y's

load('../Data/curlkick/curlkick1Y.mat')

inds{3}=[];
inds{1}=10:100;
inds{3}=find(([trials.targetcat]~=0)&([trials.disturbcat]==1));
inds{2}=find(([trials.targetcat]~=0)&([trials.disturbcat]==2));

lw=1.5;
qscale=.5;
skip=5;

for F=1:3
    for T=inds{F}
        x=100*trials(T).x;
        y=100*trials(T).y;
        forces=trials(T).f;
        
        x=[x(:,1)-center(1)+(F-1)*spacing(1), -x(:,2)+center(2)];
        plot(x(:,1),x(:,2),'color',hand)
        if F>1
            arrow(x(1:skip:end,:),x(1:skip:end,:)+qscale*forces(1:skip:end,:),gray,.3,lw);
            y=[y(:,1)-center(1)+(F-1)*spacing(1), -y(:,2)+center(2)];
            plot(y(:,1),y(:,2),'color',intent)
        end
    end
end


%% Text labels
textInvariant=['\fontname{',font,'} \fontsize{',num2str(fontsize),'} '];

%15 N cm
l0=-20;
arrow([-yw l0],[-yw l0+15*qscale],gray,.3,lw)
text(-yw,l0,[textInvariant,'15 N'],'horizontalalignment','left','verticalalignment','bottom','rotation',90)

%15 cm
l0=-20;
plot([-yw -yw+15],l0+[0 0],'k','linewidth',legendThickness)
text(-yw,l0,[textInvariant,'15 cm'],'horizontalalignment','left','verticalalignment','top')

%Experiment Labels
labs={'Baseline','Left Curl','Right Curl'};
for k=1:length(labs)
    txt=[textInvariant,labs{k}];
    text(spacing*(k-1),15,txt,'HorizontalAlignment','center','VerticalAlignment','bottom','interpreter','tex')
end

set(gca,'position',[0 0 1 1])
axis equal
%xl=xlim;
%xlim([-leftMargin xl(2)+2])

set(gcf,'color','w')
axis off