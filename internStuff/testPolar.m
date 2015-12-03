clc
clear all

nsubs=8;
exp1=2:9;
exp2=23:30;
subs=[exp1 exp2];
tlist=(1:480)';
b=floor((tlist-1)/96)+1;

%It's possible to start 2cm away from R=0
Rs=logspace(log10(2),log10(15),10);
thetabins=linspace(0,2*pi,16*3); %3 targs, 16 reaches per targ.

figure(1)
clf
hold on

largest=0;

for S=1
    if S<20
        prefix='intern';
    else
        prefix='output';
    end
    load(['./Data/',prefix,num2str(S),'.mat'])
    ind=find([trials.targcat]==0,1,'first');
    center=trials(ind).targ;
    for targs=1:3
        indt=find([trials.targcat]==targs,1,'first');
        [theta,r]=cart2pol(trials(indt).targ(1)-center(1),trials(indt).targ(2)-center(2));
        [targs theta*180/pi+180+30]
        %With this transformation, theta*180/pi+180+30, you should be able
        %to SORT by transTheta and get an ordered list of targets. BUT it
        %may make more sense to replace that 30 with the most extreme 1 as
        %the wrap point.
    end
    
    thetas=zeros(16*3*4+16*3-1,length(Rs));
    thetasy=zeros(16*3*4+16*3-1,length(Rs));
    for k=1:length(trials)
        if trials(k).targcat==0
            continue
        end
        
        %Make a trials by distances matrix of thetas for hand.
        [trials(k).theta,trials(k).r]=cart2pol(trials(k).x(:,1)-center(1),trials(k).x(:,2)-center(2));
        trials(k).r=trials(k).r*100;
        trials(k).theta=trials(k).theta*180/pi+180+30;
        for r_iter=1:length(Rs)
            [val,ind]=min(abs(trials(k).r-Rs(r_iter)));
            thetas(k,r_iter)=trials(k).theta(ind);
        end
        
        %Make a trials by distances matrix of thetas for intent.
        [trials(k).thetay,trials(k).ry]=cart2pol(trials(k).y(:,1)-center(1),trials(k).y(:,2)-center(2));
        trials(k).ry=trials(k).ry*100;
        trials(k).thetay=trials(k).thetay*180/pi+180+30;
        for r_iter=1:length(Rs)
            [val,ind]=min(abs(trials(k).ry-Rs(r_iter)));
            thetasy(k,r_iter)=trials(k).thetay(ind);
        end
        
    end
    
    targs=[trials.targcat]';
    boolNotZero=targs~=0;
    
    for block=1:5
        f=find((b==block)&boolNotZero);
        targlabs=targs(f);
        for r_iter=1:length(Rs)
            [block r_iter]
            [sorted,inds]=sort(thetas(f,r_iter));
            [targlabs(inds) sorted]
        end
    end
    
end