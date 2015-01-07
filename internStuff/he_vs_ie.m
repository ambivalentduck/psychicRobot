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
    hold on
    for kk=1:3
        plot(0,0,colors(kk))
    end
    
    for exnum=[2]
        S=subs(k,exnum);
        load(['./Data/output',num2str(S),'.mat'])
        for T=5:480
            B=floor(T/96)+1;
            if trials(T).targcat~=0
                t=trials(T).t-trials(T).t(1);
                inds=100:300;
                plot(abs(trials(T).xrot(inds,2)),abs(trials(T).yrot(inds,2)),[colors(blocksetup(B,3)),'.'],'markersize',.001);
            end
        end
    end
end

l=.05;
lo=0;
up=1;
%xlim(l*[lo up])
%ylim(l*[lo up])
plot(l*[lo up],l*[lo up],'m','linewidth',3)

legend('No forces','Forces','Forces+Intent')
axis equal
ylabel('Intent error, m')
xlabel('Hand error , m')

