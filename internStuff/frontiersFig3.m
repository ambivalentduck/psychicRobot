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
labheight=3;
asth=2.3;
toff=-2.5; %-1.5
xlift=.2;
asts=8;

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

if ~exist('stiffness.mat','file')
    for S=subs
        if S<20
            prefix='intern';
        else
            prefix='output';
        end
        load(['./Data/',prefix,num2str(S),'.mat'])
        if ~isfield(trials,'xrot')||1
            coeffs=zeros(480,4);
            for T=2:480
                p=trials(T).x;
                x=[p(:,1)-trials(T).orig(1), p(:,2)-trials(T).orig(2)];
                trials(T).start=find(vecmag(x)>.02,1,'first');
                
                trials(T).xrot=rotateProgressError(trials(T).x,trials(T).orig,trials(T).targ);
                trials(T).frot=rotateProgressError(trials(T).f,trials(T).orig,trials(T).targ);
                trials(T).yrot=rotateProgressError(trials(T).y,trials(T).orig,trials(T).targ);
                reginds=trials(T).start+(0:39);
                mgt=mean(gradient(trials(T).t(reginds)));
                
                v_error=gradient(trials(T).xrot(reginds,2))/mgt;
                a_error=gradient(v_error)/mgt;
                %trials(T).regressX=[trials(T).xrot(reginds,2) v_error a_error]; %#ok<*SAGROW>
                trials(T).regressX=[trials(T).xrot(reginds,2) v_error a_error ones(size(a_error))]; %#ok<*SAGROW>
                
                fshift=-16; %80 ms
                if reginds(1)+fshift>=1
                    trials(T).regressF=trials(T).frot(reginds+fshift,2);
                else
                    undershoot=reginds(1)+fshift;
                    trials(T).regressF=[-trials(T-1).frot(end+undershoot:end,2);trials(T).frot(1:(reginds(end)+fshift),2)];
                end
                coeffs(T,:)=trials(T).regressX\trials(T).regressF;
            end
            save(['./Data/output',num2str(S),'.mat'],'trials','params')
        end
        inds=3:2:479;
        block=floor(inds/96)+1;
        coeffs=coeffs(inds,:);
        trls=trials(inds);
        ylabels={'Stiffness, N/m','Damping, N*s/m','Mass, kg'};
        b=zeros(4,5);
        
        %figure(S)
        %clf
        for B=1:5
            state=vertcat(trls(block==B).regressX);
            force=vertcat(trls(block==B).regressF);
            [b(:,B),bint{B}]=regress(force,state);
            %[b{B},bint{B}]=regress(force,state)
            %mdl=fitlm(state,force)
            for k=1:4
                %subplot(1,4,k)
                %hold on
                %plot(B*ones(sum(block==B),1)-.4*rand(sum(block==B),1),coeffs(block==B,k),'b.')
                if k==1
                    stiff(S,B)=median(coeffs(block==B,k));
                end
                %plot(B-.2,stiff(S,B),'rx')
            end
        end
    end
    save('stiffness.mat','stiff')
else
    load('stiffness.mat')
end

stiff=abs(stiff)/100; %N/cm, easier to plot

%% Setup figure and axes
f=figure(3);
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

for E=1:2
    subplot(sph(1+E))
    hold on
    inds=Einds{E};
    rn=(jitter/8)*linspace(-1,1,8);
    
    medS=mean(stiff(inds,:));
    medS(3)-medS(2)
    plot(1:5,medS,'color',hand,'linewidth',medwidth)
    for s=1:8
        plot((1:5)-rn(s),stiff(inds(s),:),'.','markersize',msize,'color',hand)
    end
    
    xlim([.5 5.5])
    set(gca,'xtick',[])
    
    %set(gca,'xcolor','w')
    ylabel([axlabInvariant,'Stiffness (N/cm)'])
    ylim([0 2.25])
    set(gca,'ytick',0:2)
    yl=ylim;
    xl=xlim;
    text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),[textInvariant,letter{E}],'FontWeight','Bold')
    text(toff,mean(yl),mLineLabel{E},'rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','interpreter','tex')
end

%% Comparisons

subplot(sph(4))
hold on

load modelstiff.mat
ms=mean(modelstiff')';
stiff=[stiff ms/100];

k=1;
kstep=.5;
kleap=1;

xtick=[];
xticklabs={};
c=k;
for it=1:4
    xtick(end+1)=c;
    c=c+kstep;
    xtick(end+1)=c;
    c=c+kleap;
    xticklabs{end+1}='S';
    xticklabs{end+1}='R';
end

asts=8;

plot([xtick(1)-kstep xtick(end)+kstep],[0 0],'color',.7*[1 1 1])

%Hypothesis 1: I estimation innaccury
%Hypothesis 2: Performance in noise
%Hypothesis 3: Performance with I feedback
%Hypothesis 4: Which is better?
nums=[2 1;3 2; 4 3; 3 6];

for H=1:4
    pS=miniplot(stiff(exp1,nums(H,1)),stiff(exp1,nums(H,2)),k,'k.',msize,rn,asth,asts);
    ctext=['\color[rgb]{',num2str(hand(1)),',',num2str(hand(2)),',',num2str(hand(3)),'}'];
    if nums(H,2)~=6
        lab=[textInvariant,ctext,'K_',num2str(nums(H,1)),'-K_',num2str(nums(H,2))];
    else
        lab=[textInvariant,ctext,'K_',num2str(nums(H,1)),'-\color[rgb]{',num2str(intent(1)),',',num2str(intent(2)),',',num2str(intent(3)),'}K_S'];
    end
    text(k+kstep/2,labheight,lab,'HorizontalAlignment','Center')
    k=k+kstep;
    pR=miniplot(stiff(exp2,nums(H,1)),stiff(exp2,nums(H,2)),k,'k.',msize,rn,asth,asts);
    k=k+kleap;
    disp(['K',num2str(nums(H,1)),'-K',num2str(nums(H,2)),' pS=',num2str(pS),'; pR=',num2str(pR)])
end

ylabel([axlabInvariant,'Difference (N/cm)'])

set(gca,'xtick',xtick)
set(gca,'xticklabel',xticklabs)
xlim([xtick(1)-kstep xtick(end)+kstep])
ylim([-1.5 2.5])

yl=ylim;
xl=xlim;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),[textInvariant,'D'],'FontWeight','Bold')