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
            T
            B=floor(T/96)+1;
            if trials(T).targcat~=0
                if B==5
                    continue
                end
                t=trials(T).t-trials(T).t(1);
                inds=1:min(350,length(t));
                shand=sign(trials(T).yrot(inds,2));
                plot(t(inds),log10(abs(trials(T).xrot(inds,2))./abs(trials(T).yrot(inds,2))),[colors(blocksetup(B,3)),'.'],'markersize',.001);
            end
        end
    end
end
%axis equal
%legend('No forces','Forces','Forces+Intent')
ylabel('log(F/a)')
xlabel('time, s')

