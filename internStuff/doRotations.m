clc
clear all
close all

subs=[2:8 23:29];

for S=subs
    if S<20
        prefix='intern';
    else
        prefix='output';
    end
    load(['./Data/',prefix,num2str(S),'.mat'])
    if ~isfield(trials,'xrot')|1
        for T=2:480
            if (T>(2*96))&&(T<=(96*3))
                p=trials(T).y;
            else
                p=trials(T).x;
            end
            x=[p(:,1)-trials(T).orig(1), p(:,2)-trials(T).orig(2)];
            trials(T).start=max(find(vecmag(x)>.02,1,'first'),1);
            
            trials(T).xrot=rotateProgressError(trials(T).x,trials(T).orig,trials(T).targ);
            trials(T).frot=rotateProgressError(trials(T).f,trials(T).orig,trials(T).targ);
            trials(T).yrot=rotateProgressError(trials(T).y,trials(T).orig,trials(T).targ);
        end
        save(['./Data/output',num2str(S),'.mat'],'trials','params')
    end
end