clc
clear all

%% Set constants and initialize

%Measurements in cm. Movements are 15 cm long

[hand,intent]=getColorScheme;

H=20;
W=18; %Do NOT adjust

yw=15*cosd(30);

leftMargin=60;
legendThickness=2;
legendSpacing=7;
legendp=[-yw 61];
center=[0 48];
spacing=[40 30];
font='Arial';
fontsize=12;
interpreter='tex';
handfudge=0;
intentfudge=0;

snums=[4 28]; %Exp1=2-9, Exp2=23-30; 4 and 28 are the easiest visually

f=figure(1);
clf
set(f,'units','centimeters');
set(f,'Position',[3,3,18,8])

hold on

tlist=(1:480)';
b=floor((tlist-1)/96)+1;

%% Draw Y's
for E=1:2
    load(['./Data/output',num2str(snums(E)),'.mat'])
    for T=3:480
        if trials(T).targcat==0
            continue
        end
        B=b(T);
        x=100*trials(T).x;
        y=100*trials(T).y;
        
        x=[x(:,1)-center(1)+(B-1)*spacing(1), -x(:,2)+center(2)+(2-E)*spacing(2)];
        plot(x(:,1),x(:,2),'color',hand)
        if B==3
            y=[y(:,1)-center(1)+(B-1)*spacing(1), -y(:,2)+center(2)+(2-E)*spacing(2)];
            plot(y(:,1),y(:,2),'color',intent)
        end
        
    end
end

%% Text labels

textInvariant=['\fontname{',font,'} \fontsize{',num2str(fontsize),'} '];
%Forces
txt=[textInvariant,'Forces On'];
text(legendp(1),legendp(2),txt,'HorizontalAlignment','right','VerticalAlignment','middle','interpreter','tex')
plot([-yw+spacing(1)*1,yw+spacing(1)*3],legendp(2)+[0 0],'k','linewidth',legendThickness)

%Cursor=hand
txt=[textInvariant,'Cursor = \color[rgb]{',num2str(hand(1)),',',num2str(hand(2)),',',num2str(hand(3)),'}Hand'];
text(legendp(1)+handfudge,legendp(2)-legendSpacing,txt,'HorizontalAlignment','right','VerticalAlignment','middle','interpreter','tex')
plot([-yw yw+spacing(1)*1],legendp(2)-legendSpacing+[0 0],'k','linewidth',legendThickness)
plot([-yw+spacing(1)*3 yw+spacing(1)*4],legendp(2)-legendSpacing+[0 0],'k','linewidth',legendThickness)

%Cursor=intent
txt=[textInvariant,'Cursor = \color[rgb]{',num2str(intent(1)),',',num2str(intent(2)),',',num2str(intent(3)),'}Intent'];
text(legendp(1)+intentfudge,legendp(2)-2*legendSpacing,txt,'HorizontalAlignment','right','VerticalAlignment','middle','interpreter','tex')
plot([-yw yw]+spacing(1)*2,legendp(2)-2*legendSpacing+[0 0],'k','linewidth',legendThickness)

%15 cm
l0=-20;
plot([-yw -yw+15],l0+[0 0],'k','linewidth',legendThickness)
text(-yw+7.5,l0,[textInvariant,'15 cm'],'horizontalalignment','center','verticalalignment','top')

%Experiment Labels
yoff=15*(sind(30)-1)/2;
txt={[textInvariant,'Standard']};
text(legendp(1),spacing(2)+yoff,txt,'HorizontalAlignment','right','VerticalAlignment','middle','interpreter','tex')

txt={[textInvariant,'Reduced'],[textInvariant, 'Stiffness']};
text(legendp(1),yoff,txt,'HorizontalAlignment','right','VerticalAlignment','middle','interpreter','tex')


set(gca,'position',[0 0 1 1])
axis equal
xl=xlim;
xlim([-leftMargin xl(2)+2])

set(gcf,'color','w')
axis off