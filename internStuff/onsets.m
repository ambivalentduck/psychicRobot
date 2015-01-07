clc
clear all
close all

subs=[2:8; 23:29]';
blocksetup=[0 0 1;1 0 2; 1 1 3; 1 0 2; 0 0 1];
colors='krg';
markers={'-','-.','--'};

for k=2 %:7
    figure(k)
    clf
    subplot(3,1,1)
    hold on
    subplot(3,1,2)
    hold on
    subplot(3,1,3)
    hold on
    for kk=1:3
        plot(0,0,colors(kk))
    end
    
    for exnum=[2]
        S=subs(k,exnum);
        load(['./Data/output',num2str(S),'.mat'])
        for T=5:480
            T
            B=floor(T/96)+1;
            if B==5
                continue
            end
            if trials(T).targcat~=0
                t=trials(T).t-trials(T).t(1);
                inds=1:min(350,length(t));
                shand=sign(trials(T).yrot(inds,2));
                subplot(3,1,1)
                plot(t(inds),vecmag(trials(T).v(inds,:)),[colors(blocksetup(B,3)),'-'],'markersize',.001);
                subplot(3,1,2)
                plot(t(inds),vecmag(trials(T).f(inds,:)),[colors(blocksetup(B,3)),'-'],'markersize',.001);
                subplot(3,1,3)
                plot(t(inds),log10(abs(trials(T).xrot(inds,2))./abs(trials(T).yrot(inds,2))),[colors(blocksetup(B,3)),'.'],'markersize',.001)
            end
        end
    end
    plot([0 1.5],[0 0],'m')
end
%axis equal
%legend('No forces','Forces','Forces+Intent')
subplot(3,1,2)
ylabel('force, N')
xlabel('time, s')

subplot(3,1,1)
ylabel('speed, m/s')

