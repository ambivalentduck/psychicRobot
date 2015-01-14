clc
clear all

t=0:.01:1;

figure(1)
clf

outsum=zeros(length(t),3);

tc1=.3
tc=[tc1 .5 2*.5-tc1];
ts=[(.5-tc(1))*2 .4 (tc(3)-.5)*2];

colors='rgb';

for k=1:3
    [kern{k},out{k}]=slmj5op(t,tc(k),ts(k));
    if k~=2
        out{k}=out{k}/1;
    end
    outsum=outsum+out{k};
    
    for p=1:3
        subplot(3,1,p)
        hold on
        plot(t,out{k}(:,p),colors(k))
    end
end

for p=1:3
    subplot(3,1,p)
    plot(t,outsum(:,p),'k')
end

tc1s=0:.01:.5;
sums=zeros(length(t),length(tc1s));
for c=1:length(tc1s)
    tc1=tc1s(c);
    tc=[tc1 .5 2*.5-tc1];
    ts=[(.5-tc(1))*2 .4 (tc(3)-.5)*2];
    
    for k=1:3
        kern=slmj5op(t,tc(k),ts(k));
        if k~=2
            kern=kern/10;
        end
        sums(:,c)=sums(:,c)+kern;
    end
end