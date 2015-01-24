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
    b=floor((tlist-1)/96)+1
    t=0:.005:149*.005;
    for treat=1:3
        f=find((tlist>5)&mod(tlist,2)&(blocksetup(b,3)==treat));
        hc=horzcat(tr(f).qoi);
        mhc{k,treat}=mean(hc');
        plot(t,mhc{k,treat},'-','color',.5*colors(treat,:),'linewidth',5)
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
hold on

tls=5:2:479
for k=1:7
    plot(
end
legend('No forces','Forces','Forces+Intent')
ylabel('|Hand Error|-|Intent Error|, m')
xlabel('time, s')
plot([0 t(end)],[0 0],'k-')
xlim([0 t(end)])

