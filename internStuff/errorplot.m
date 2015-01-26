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

fh=figure(18);
clf

subplot(3,4,[2:4 6:8])
hold on
exp1=2:9;
exp2=23:30;
rn=.1*(rand(nsubs,1)-.5);

spacer=5;

for s=1:nsubs
    plot((1:5)+rn(s),jerror1(exp1(s),:)+spacer,'r-')
    plot((1:5)+rn(s),jerror1(exp1(s),:)+spacer,'r.')
    plot((1:5)+rn(s),jerror1(exp2(s),:)+spacer,'b-')
    plot((1:5)+rn(s),jerror1(exp2(s),:)+spacer,'b.')
    
    plot((1:5)+rn(s),jerror2(exp1(s),:),'r-')
    plot((1:5)+rn(s),jerror2(exp1(s),:),'r.')
    plot((1:5)+rn(s),jerror2(exp2(s),:),'b-')
    plot((1:5)+rn(s),jerror2(exp2(s),:),'b.')
end
set(gca,'xtick',1:5)
ticks=[0:3 spacer+(0:3)];
extras=spacer+4+(1:3)*.9;
set(gca,'ytick',[ticks extras])
ylabs={};
nums=[0:3 0:3];
for k=1:length(nums)
    ylabs{k}=num2str(nums(k));
end
ylabs{end+1}='Cursor=Intent';
ylabs{end+1}='Cursor=Hand';
ylabs{end+1}='Forces On';
set(gca,'yticklabel',ylabs)

xlim([.5 5.5])
uyl=ylabel('Hand');
g=get(uyl,'Position');
set(uyl,'Position',[.25 g(2:3)])
tap=[.1277,.77,0,0];
ann=annotation(fh,'textarrow',[.5 .5],[.5 .5],'String','Max Deviation (cm)','TextRotation',90,'HeadStyle','none','LineStyle', 'none');
set(ann,'Position',tap)
ann2=annotation(fh,'textarrow',[.5 .5],[.5 .5],'String','Intent','TextRotation',90,'HeadStyle','none','LineStyle', 'none');
tap2=[.1277,.75,0,0];
set(ann2,'Position',tap2)
plot([1.5 4.5],extras(3)*[1 1],'k')
plot([.5 2.5],extras(2)*[1 1],'k',[3.5 5.5],extras(2)*[1 1],'k')
plot([2.5 3.5],extras(1)*[1 1],'k')
ylim([0 extras(3)+.25])

xlabel('Block')


xl=[.7 2.3];
%Hypothesis 1: I estimation innaccury
msize=10;
subplot(3,4,9)
miniplot(jerror1(exp1,1),jerror2(exp1,1),jerror1(exp2,1),jerror2(exp2,1),msize,xl,rn)
ylabel('\Delta Max Deviation (cm)')
title('H_1-I_1')

%Hypothesis 2: Performance in noise
subplot(3,4,10)
miniplot(jerror1(exp1,2),jerror2(exp1,2),jerror1(exp2,2),jerror2(exp2,2),msize,xl,rn)
title('H_2-I_2')

%Hypothesis 3: Performance with I feedback
subplot(3,4,11)
miniplot(jerror1(exp1,3),jerror2(exp1,3),jerror1(exp2,3),jerror2(exp2,3),msize,xl,rn)
title('H_3-I_3')

%Hypothesis 4: Which is better?
subplot(3,4,12)
miniplot(jerror1(exp1,2),jerror2(exp1,3),jerror1(exp2,2),jerror2(exp2,3),msize,xl,rn)
title('H_2-I_3')
xlabel('Experiment')
