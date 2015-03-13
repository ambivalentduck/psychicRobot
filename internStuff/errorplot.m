clc
clear all

nsubs=8;
subs=[2:9 23:30];
tlist=(1:480)';
b=floor((tlist-1)/96)+1;

if ~exist('jerror.mat','file')
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
            if S==5
            jerror1(S,blocks)=mean(v(v(:,1)<.05,1));
            jerror2(S,blocks)=mean(v(v(:,2)<.05,2));    
            else
            m=mean(v);
            jerror1(S,blocks)=m(1);
            jerror2(S,blocks)=m(2);
            end
        end
    end
    save('jerror.mat','jerror1','jerror2')
else
    load('jerror.mat')
end

jerror1=jerror1*100;
jerror2=jerror2*100;

%% New stuff
msize=8;

f=figure(3);
set(f,'units','in');
width=6;
lmargin=1;
rmargin=.05;

set(f,'Position',[5,5,3,5.5])
clf

margin=.1;
ploth=1;

sp=subplot('Position',[0 .9 1 .05]);
set(sp,'Units','in','Position',[lmargin 4.75 3-rmargin-lmargin .5])
hold on
plot([1.5 4.5],3*[1 1],'k')
plot([2.5 3.5],1*[1 1],'k')
plot([.5 2.5],2*[1 1],'k',[3.5 5.5],2*[1 1],'k')

set(gca,'xtick',[])
set(gca,'ytick',1:3)
ylabs={'Cursor=Intent','Cursor=Hand','Forces On'};

set(gca,'yticklabel',ylabs)
set(gca,'xcolor','w')

ylim([.5 3.5])

yl=ylim;
xl=xlim;
xrat=-.9/1.95;
yrat=1+.1/.5;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),'A','FontWeight','Bold')

%% Hand

sp1=subplot('Position',[0 .4 1 .2])
set(sp1,'Units','in','Position',[lmargin 3.5 3-rmargin-lmargin 1])
hold on
exp1=2:9;
exp2=23:30;
rn=.1*linspace(-1,1,8);

for s=1:8
    plot((1:5)-rn(s),jerror1(exp1(s),:),'r-')
    plot((1:5)-rn(s),jerror1(exp1(s),:),'r.')
    plot((1:5)-rn(s),jerror1(exp2(s),:),'b-')
    plot((1:5)-rn(s),jerror1(exp2(s),:),'b.')
end
xlim([.5 5.5])
set(gca,'xtick',[])
%set(gca,'xcolor','w')
ylabel('Hand Error (cm)')
yl=ylim;
xl=xlim;
xrat=-.9/1.95;
yrat=1+.1/1;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),'B','FontWeight','Bold')

%% Intent
sp1=subplot('Position',[0 .2 1 .3])
set(sp1,'Units','in','Position',[lmargin 2.25 3-rmargin-lmargin 1])
hold on
exp1=2:9;
exp2=23:30;
rn=.1*linspace(-1,1,8);


for s=1:8
    plot((1:5)-rn(s),jerror2(exp1(s),:),'r-')
    plot((1:5)-rn(s),jerror2(exp1(s),:),'r.')
    plot((1:5)-rn(s),jerror2(exp2(s),:),'b-')
    plot((1:5)-rn(s),jerror2(exp2(s),:),'b.')
end

xlim([.5 5.5])
ylabel('Intent Error (cm)')
xlabel('Block')
yl=ylim;
xl=xlim;
xrat=-.9/1.95;
yrat=1+.1/1;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),'C','FontWeight','Bold')

%% Comparisons

sp=subplot('Position',[0 0 1 .1])
set(sp,'Units','in','Position',[lmargin .5 3-rmargin-lmargin 1])

hold on

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
    xticklabs{end+1}='1';
    xticklabs{end+1}='2';
end

asts=8;
asth=1.5;

plot([xtick(1)-kstep xtick(end)+kstep],[0 0],'color',.7*[1 1 1])

%Hypothesis 1: forces increase stiffness (2-1)

labheight=2.2;

%Hypothesis 1: I estimation innaccury
miniplot(jerror1(exp1,1),jerror2(exp1,1),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'H_1-I_1','HorizontalAlignment','Center')
k=k+kstep;
miniplot(jerror1(exp2,1),jerror2(exp2,1),k,'r.',msize,rn,asth,asts)
k=k+kleap;

%Hypothesis 2: Performance in noise
miniplot(jerror1(exp1,2),jerror2(exp1,2),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'H_2-I_2','HorizontalAlignment','Center')
k=k+kstep;
miniplot(jerror1(exp2,2),jerror2(exp2,2),k,'r.',msize,rn,asth,asts)
k=k+kleap;

%Hypothesis 3: Performance with I feedback
miniplot(jerror1(exp1,3),jerror2(exp1,3),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'H_3-I_3','HorizontalAlignment','Center')
k=k+kstep;
miniplot(jerror1(exp2,3),jerror2(exp2,3),k,'r.',msize,rn,asth,asts)
k=k+kleap;

%Hypothesis 4: Which is better?
miniplot(jerror1(exp1,2),jerror2(exp1,3),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'H_2-I_3','HorizontalAlignment','Center')
k=k+kstep;
miniplot(jerror1(exp2,2),jerror2(exp2,3),k,'r.',msize,rn,asth,asts)
k=k+kleap;

ylabel('\Delta Max Deviation (cm)')

set(gca,'xtick',xtick)
set(gca,'xticklabel',xticklabs)
xlim([xtick(1)-kstep xtick(end)+kstep])

xlabel('Experiment')

yl=ylim;
xl=xlim;
xrat=-.9/1.95;
yrat=1+.1/1;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),'D','FontWeight','Bold')


