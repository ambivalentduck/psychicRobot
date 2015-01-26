clc
clear all
close all

nsubs=8;
subs=[2:9 23:30];
tlist=(1:480)';
b=floor((tlist-1)/96)+1;


for S=5 %subs
    prefix='output';
    
    load(['./Data/',prefix,num2str(S),'.mat'])
    
    for T=5:480
        if trials(T).targcat~=0
            inds=trials(T).start:min(trials(T).start+49,length(trials(T).t));
            tr(T).jerrs=max([abs(trials(T).xrot(inds,2)) abs(trials(T).yrot(inds,2))]);
            if any(tr(T).jerrs>.07)
                figure(T)
                clf
                hold on
                plot(trials(T).xrot(:,1),trials(T).xrot(:,2),'b')
                plot(trials(T).yrot(:,1),trials(T).yrot(:,2),'r')
                axis equal
            end
        end
    end

end
