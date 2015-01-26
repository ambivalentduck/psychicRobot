clc
clear all
close all

subs=[2:8; 23:29]';
blocksetup=[0 0 1;1 0 2; 1 1 3; 1 0 2; 0 0 1];
colors=[0 0 1;1 0 0; 0 1 0];
markers={'-','-.','--'};

for k=1:7
    figure(k)
    clf
    subplot(3,1,3)
    hold on
    for kk=1:3
        plot(0,0,'color',colors(kk,:))
    end
    
    for exnum=[2]
        S=subs(k,exnum);
        load(['./Data/output',num2str(S),'.mat'])
        tlength=500*ones(480,1);
        for T=5:480
            B=floor(T/96)+1;
            
            if trials(T).targcat~=0
                inds=trials(T).start:min(trials(T).start+400,length(trials(T).t));
                t=trials(T).t(inds)-trials(T).t(trials(T).start);
                tlength(T)=length(inds);
                qoi=abs(trials(T).xrot(inds,2))-abs(trials(T).yrot(inds,2));
                subplot(3,1,3)
                plot(t,qoi,'-','color',colors(blocksetup(B,3),:),'markersize',.001);
                subplot(3,1,2)
                hold on
                plot(t,vecmag(trials(T).f(inds,:)),'-','color',colors(blocksetup(B,3),:))
                subplot(3,1,1)
                hold on
                plot(t,vecmag(trials(T).v(inds,:)),'-','color',colors(blocksetup(B,3),:))
                tr(T).qoi=zeros(150,1);
                lqoi=min(length(qoi),150);
                tr(T).qoi(1:lqoi)=qoi(1:lqoi);
                tr(T).jerrs=mean([abs(trials(T).xrot(inds(1:50),2)) abs(trials(T).yrot(inds(1:50),2))]);
                rmsx(k,T)=sqrt(mean(trials(T).xrot(inds,2).^2));
                rmsy(k,T)=sqrt(mean(trials(T).yrot(inds,2).^2));
            end
        end
    end
    %axis equal
    subplot(3,1,3)
    legend('No forces','Forces','Forces+Intent')
    ylabel('|Hand Error|-|Intent Error|, m')
    xlabel('time, s')
    
    subplot(3,1,2)
    ylabel('|Force|')
    subplot(3,1,1)
    ylabel('Speed')
    
    subplot(3,1,3)
    tlist=(1:480)';
    b=floor((tlist-1)/96)+1;
    t=0:.005:149*.005;
    for treat=1:3
        f=find((tlist>5)&mod(tlist,2)&(blocksetup(b,3)==treat));
        hc=horzcat(tr(f).qoi);
        mhc{k,treat}=mean(hc');
        plot(t,mhc{k,treat},'-','color',.5*colors(treat,:),'linewidth',5)
    end
    
    for blocks=1:5
        f=find((tlist>5)&mod(tlist,2)&(b==blocks));
        jerror{k,blocks}=mean(vertcat(tr(f).jerrs));
    end
end

figure(17)
clf
hold on

for k=1:7
    for treat=1:3
        plot(t,mhc{k,treat},'-','color',colors(treat,:),'linewidth',1)
    end
end
legend('No forces','Forces','Forces+Intent')
ylabel('|Hand Error|-|Intent Error|, m')
xlabel('time, s')
plot([0 t(end)],[0 0],'k-')
xlim([0 t(end)])

figure(18)
clf

subplot(2,4,2:4)
hold on
exp1=2:8;
exp2=23:29;
rn=.1*(rand(7,1)-.5);


for s=1:7
    plot((1:5)-rn(s),jerror(exp1(s),:),'r-')
    plot((1:5)-rn(s),jerror(exp1(s),:),'r.')
    plot((1:5)-rn(s),jerror(exp2(s),:),'b-')
    plot((1:5)-rn(s),jerror(exp2(s),:),'b.')
end
set(gca,'xtick',1:5)
set(gca,'ytick',[0:50:200 225 250 275])
ylabs={};
for k=0:50:200
    ylabs{k/50+1}=num2str(k);
end
ylabs{end+1}='Cursor=Deduced Intent';
ylabs{end+1}='Cursor=Hand';
ylabs{end+1}='Forces On';
length(ylabs)
length([0:50:200 225 250 275])
set(gca,'yticklabels',ylabs)
ylim([0 300])

xlim([.5 5.5])
ylabel('jerrorness, N/m')
xlabel('Block')
plot([1.5 4.5],275*[1 1],'k')
plot([2.5 3.5],225*[1 1],'k')
plot([.5 2.5],250*[1 1],'k',[3.5 5.5],250*[1 1],'k')

xl=[.7 2.3];
%Hypothesis 1: forces increase jerrorness (2-1)
msize=3;
subplot(2,4,5)
hold on
plot(1+rn,jerror(exp1,2)-jerror(exp1,1),'r.','markersize',msize)
plot(2+rn,jerror(exp2,2)-jerror(exp2,1),'b.','markersize',msize)
set(gca,'xtick',[1 2])
xlim(xl)
ylabel('\Deltajerrorness, N/m')
title('2-1')

%Hypothesis 2: intent decreases jerrorness despite forces (3-2)
subplot(2,4,6)
hold on
plot(1+rn,jerror(exp1,3)-jerror(exp1,2),'r.','markersize',msize)
plot(2+rn,jerror(exp2,3)-jerror(exp2,2),'b.','markersize',msize)
set(gca,'xtick',[1 2])
xlim(xl)
title('3-2')

%Hypothesis 2.5: indistinguishable from baseline (3-1)
subplot(2,4,7)
hold on
plot(1+rn,jerror(exp1,3)-jerror(exp1,1),'r.','markersize',msize)
plot(2+rn,jerror(exp2,3)-jerror(exp2,1),'b.','markersize',msize)
set(gca,'xtick',[1 2])
xlim(xl)
title('3-1')
%Hypothesis 3: jerrorness increases again when intent is removed (4-3)
subplot(2,4,8)
hold on
plot(1+rn,jerror(exp1,4)-jerror(exp1,3),'r.','markersize',msize)
plot(2+rn,jerror(exp2,4)-jerror(exp2,3),'b.','markersize',msize)
set(gca,'xtick',[1 2])
xlim(xl)
title('4-3')

