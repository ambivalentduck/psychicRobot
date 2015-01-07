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
    subplot(2,1,1)
    hold on
    for kk=1:3
        plot(0,0,colors(kk))
    end
    subplot(2,1,2)
    hold on
    
    for exnum=[2]
        S=subs(k,exnum);
        load(['./Data/output',num2str(S),'.mat'])
        try
            error
            trials(3).xrot+2;
        catch
            for T=1:480
                trials(T).arot=rotateProgressError(trials(T).a,trials(T).x(100,:),trials(T).targ);
                trials(T).xrot=rotateProgressError(trials(T).x,trials(T).x(100,:),trials(T).targ);
                trials(T).frot=rotateProgressError(trials(T).f,trials(T).x(100,:),trials(T).targ);
                trials(T).yrot=rotateProgressError(trials(T).y,trials(T).y(100,:),trials(T).targ);
            end
            save(['./Data/output',num2str(S),'.mat'],'trials','params')
        end
        for T=5:480
            B=floor(T/96)+1;
            if trials(T).targcat~=0
                t=trials(T).t-trials(T).t(1);
                subplot(2,1,1)
                plot(t,abs(trials(T).frot(:,2))./abs(trials(T).xrot(:,2)),colors(blocksetup(B,3)),'linewidth',.001);
                subplot(2,1,2)
                plot(t,abs(trials(T).frot(:,2))./abs(trials(T).yrot(:,2)),colors(blocksetup(B,3)),'linewidth',.001);
            end
        end
    end
end

subplot(2,1,1)
legend('No forces','Forces','Forces+Intent')
title('Hand')
xlim([0 2])
subplot(2,1,2)
title('Extraction')
xlim([0 2])
