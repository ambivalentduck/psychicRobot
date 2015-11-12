clc
clear all

%% Set constants
font='Arial';
fontsize=12;
axlabsize=10;
interpreter='tex';
textInvariant=['\fontname{',font,'} \fontsize{',num2str(fontsize),'} '];
axlabInvariant=['\fontname{',font,'} \fontsize{',num2str(axlabsize),'} '];
[hand,intent]=getColorScheme;

lmargin=3;
rmargin=.25;
msize=8;
figW=8.5;
axW=8.5-lmargin-rmargin;
hiStagger=.2;
jitter=.8;
medwidth=2;
labheight=2.35;
asth=1.85;
toff=-2.5; %-1.5
xlift=.2;
asts=8;

nums=[3 3; 2 3; 4 3];

xrat=-3.2/axW;
yrat=1+.003/1;

exp1=2:9;
exp2=23:30;

spy=[13.75, 15.75;
    9.25, 13.25;
    4.5, 8.5;
    .5, 3.5];

figH=spy(1,2)+1;

%% Load data
nsubs=8;
subs=[exp1 exp2];
tlist=(1:480)';
b=floor((tlist-1)/96)+1;

p43=[];
if ~exist('jerror.mat','file')
    flatten=@mean;
    for S=subs
        prefix='output';
        
        load(['./Data/',prefix,num2str(S),'.mat'])
        
        for T=5:480
            if trials(T).targcat~=0
                inds=trials(T).start:min(trials(T).start+49,length(trials(T).t));
                tr(T).jerrs=max([abs(trials(T).xrot(inds,2)) abs(trials(T).yrot(inds,2))]);
            end
        end
        
        for blocks=1:5
            f=find((tlist>5)&mod(tlist,2)&(b==blocks));
            v=vertcat(tr(f).jerrs);
            if S==5 %One subject made beautiful reaches...to the wrong target
                yup=(v(:,1)<.05)&(v(:,2)<.05);
                v=v(yup,:); %5cm error is the magic cutoff to get wrong targets
            end
            comparison{blocks}=v;
            m=flatten(v);
            jerror1(S,blocks)=m(1);
            jerror2(S,blocks)=m(2);            
        end
        [h,p43(S)]=ttest2(comparison{3}(:,2),comparison{4}(:,1));
    end
    save('jerror.mat','jerror1','jerror2')
else
    load('jerror.mat')
end

jerror1=jerror1*100;
jerror2=jerror2*100;

%% Setup figure and axes
f=figure(2);
clf
set(f,'units','centimeters','Position',[5,5,figW,figH]);

sph=zeros(4,1);
for k=1:4
    sph(k)=subplot(4,1,5-k);
end
for k=1:4
    set(sph(k),'Units','centimeters','Position',[lmargin spy(k,1) figW-rmargin-lmargin spy(k,2)-spy(k,1)]);
end

%% Plot methods figure
subplot(sph(1))
hold on

plot([1.5 4.5],3*[1 1],'k')
plot([2.5 3.5],1*[1 1],'k')
plot([.5 2.5],2*[1 1],'k',[3.5 5.5],2*[1 1],'k')

xlim([.5 5.5])
ylim([.5 3.5])

set(gca,'xtick',1:5)
set(gca,'xaxislocation','top')
set(gca,'ytick',1:3)
set(gca,'yticklabel',{[],[],[]})

ylabs={[textInvariant,'Cursor = \color[rgb]{',num2str(intent(1)),',',num2str(intent(2)),',',num2str(intent(3)),'}Intent'],...
    [textInvariant,'Cursor = \color[rgb]{',num2str(hand(1)),',',num2str(hand(2)),',',num2str(hand(3)),'}Hand'],...
    [textInvariant,'Forces On']};
for k=1:3
    text(.4,k,ylabs{k},'interpreter','tex','HorizontalAlignment','right','VerticalAlignment','middle')
end

title([textInvariant,'Block'])
yl=ylim;
xl=xlim;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),[textInvariant,'A'],'FontWeight','Bold')

%% Experiments 1 and 2
Einds={exp1,exp2};
letter={'B','C'};
mLineLabel={{[textInvariant,'\bfS\rmtandard']},...
    {[textInvariant,'\bfR\rmeduced'],[textInvariant,'Stiffness']}};
stagger=[0 0 hiStagger/2 0 0];
for E=1:2
    subplot(sph(1+E))
    hold on
    inds=Einds{E};
    rn=(jitter/8)*linspace(-1,1,8);
    
    medX=mean(jerror1(inds,:));
    medY=mean(jerror2(inds,:));
    medC=medX;
    medC(3)=medY(3);
    
    blueline1=[1 2 2.5];
    blueline2=[3.5 4 5];
    redline=[2.5 3 3.5];
    plot(blueline1,interp1(1:5,medC,blueline1),'color',hand,'linewidth',medwidth)
    plot(blueline2,interp1(1:5,medC,blueline2),'color',hand,'linewidth',medwidth)
    plot(redline,interp1(1:5,medC,redline),'color',intent,'linewidth',medwidth)
    
    for s=1:8
        plot((1:5)-stagger-rn(s),jerror1(inds(s),:),'.','markersize',msize,'color',hand)
        plot(3+hiStagger/2-rn(s),jerror2(inds(s),3),'.','markersize',msize,'color',intent)
    end
    
    xlim([.5 5.5])
    set(gca,'xtick',[])
    
    %set(gca,'xcolor','w')
    ylabel([axlabInvariant,'Error (cm)'])
    yl=ylim;
    xl=xlim;
    text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),[textInvariant,letter{E}],'FontWeight','Bold')
    text(toff,mean(yl),mLineLabel{E},'rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','interpreter','tex')
end

%% Comparisons

subplot(sph(4))
hold on

k=1;
kstep=.5;
kleap=1;

xtick=[];
xticklabs={};
c=k;
for it=1:size(nums,1)
    xtick(end+1)=c;
    c=c+kstep;
    xtick(end+1)=c;
    c=c+kleap;
    xticklabs{end+1}='S';
    xticklabs{end+1}='R';
end

plot([xtick(1)-kstep xtick(end)+kstep],[0 0],'color',.7*[1 1 1])

%Hypothesis 1: I estimation innaccury
%Hypothesis 2: Performance in noise
%Hypothesis 3: Performance with I feedback
%Hypothesis 4: Which is better?

for H=1:size(nums,1)
    [pS,statsS,meanS]=miniplot(jerror1(exp1,nums(H,1)),jerror2(exp1,nums(H,2)),k,'k.',msize,rn,asth,asts);
    lab=[textInvariant,'\color[rgb]{',num2str(hand(1)),',',num2str(hand(2)),',',num2str(hand(3)),...
        '}H_',num2str(nums(H,1)),'\color{black}-\color[rgb]{',num2str(intent(1)),',',num2str(intent(2)),',',num2str(intent(3)),...
        '}I_',num2str(nums(H,2))];
    text(k+kstep/2,labheight,lab,'HorizontalAlignment','Center')
    k=k+kstep;
    [pR,statsR,meanR]=miniplot(jerror1(exp2,nums(H,1)),jerror2(exp2,nums(H,2)),k,'k.',msize,rn,asth,asts);
    k=k+kleap;
    disp(['H',num2str(nums(H,1)),'-I',num2str(nums(H,2))])
    disp(['pS=',num2str(pS),'; TS=',num2str(statsS.tstat),'; dfS=',num2str(statsS.df),'; meanS=',num2str(meanS),'; semS=',num2str(statsS.sd/sqrt(statsS.df))])
    disp(['pR=',num2str(pR),'; TR=',num2str(statsR.tstat),'; dfR=',num2str(statsR.df),'; meanR=',num2str(meanR),'; semR=',num2str(statsR.sd/sqrt(statsR.df))])
end

ylabel([axlabInvariant,'Difference (cm)'])

set(gca,'xtick',xtick)
set(gca,'xticklabel',xticklabs)
xlim([xtick(1)-kstep xtick(end)+kstep])

yl=ylim;
xl=xlim;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),[textInvariant,'D'],'FontWeight','Bold')


