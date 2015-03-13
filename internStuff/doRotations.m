clc
clear all


subs=[2:9 23:30];

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
%                             if (T>(2*96))&&(T<=(96*3))
%                                 p=trials(T).y;
%                             else
%                                 p=trials(T).x;
%                             end
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

%% Post-hoc analysis


msize=8;

f=figure(3);
set(f,'units','in');
width=6;
lmargin=1;
rmargin=.05;

set(f,'Position',[5,5,3,4.25])
clf

margin=.1;
ploth=1;

sp=subplot('Position',[0 .8 1 .1]);
set(sp,'Units','in','Position',[lmargin 3.5 3-rmargin-lmargin .5])
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


sp1=subplot('Position',[0 .4 1 .2])
set(sp1,'Units','in','Position',[lmargin 2.25 3-rmargin-lmargin 1])
hold on
exp1=2:9;
exp2=23:30;
rn=.1*linspace(-1,1,8);

stiff=abs(stiff);

for s=1:8
    plot((1:5)-rn(s),stiff(exp1(s),:),'r-')
    plot((1:5)-rn(s),stiff(exp1(s),:),'r.')
    plot((1:5)-rn(s),stiff(exp2(s),:),'b-')
    plot((1:5)-rn(s),stiff(exp2(s),:),'b.')
end


xlim([.5 5.5])
ylabel('Stiffness (N/m)')
xlabel('Block')
yl=ylim;
xl=xlim;
xrat=-.9/1.95;
yrat=1+.1/1;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),'B','FontWeight','Bold')

sp=subplot('Position',[0 0 1 .2])
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
asth=199;

plot([xtick(1)-kstep xtick(end)+kstep],[0 0],'color',.7*[1 1 1])

%Hypothesis 1: forces increase stiffness (2-1)

labheight=250;

miniplot(stiff(exp1,2),stiff(exp1,1),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'2-1','HorizontalAlignment','Center')
k=k+kstep;
miniplot(stiff(exp2,2),stiff(exp2,1),k,'r.',msize,rn,asth,asts)
k=k+kleap;

%Hypothesis 2: intent decreases stiffness despite forces (3-2)
miniplot(stiff(exp1,3),stiff(exp1,2),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'3-2','HorizontalAlignment','Center')
k=k+kstep;
miniplot(stiff(exp2,3),stiff(exp2,2),k,'r.',msize,rn,asth,asts)
k=k+kleap;

%Hypothesis 3: stiffness increases again when intent is removed (4-3)
miniplot(stiff(exp1,4),stiff(exp1,3),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'4-3','HorizontalAlignment','Center')
k=k+kstep;
miniplot(stiff(exp2,4),stiff(exp2,3),k,'r.',msize,rn,asth,asts)
k=k+kleap;

%Hypothesis 4: people are essentially adopting model stiffness
load modelstiff.mat
ms=mean(modelstiff')';
miniplot(stiff(exp1,3),ms(exp1),k,'b.',msize,rn,asth,asts)
text(k+kstep/2,labheight,'3-Model','HorizontalAlignment','Center')
k=k+kstep;
miniplot(stiff(exp2,3),ms(exp2),k,'r.',msize,rn,asth,asts)
k=k+kleap;


ylabel('\DeltaStiffness (N/m)')

set(gca,'xtick',xtick)
set(gca,'xticklabel',xticklabs)
xlim([xtick(1)-kstep xtick(end)+kstep])

xlabel('Experiment')

yl=ylim;
xl=xlim;
xrat=-.9/1.95;
yrat=1+.1/1;
text(xrat*xl(2)+(1-xrat)*xl(1),yrat*yl(2)+(1-yrat)*yl(1),'C','FontWeight','Bold')
