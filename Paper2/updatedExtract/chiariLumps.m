clc
clear all
load('pulse1.mat');
close all

N=14;
onset=find(vecmag(trials(N).v)>.05,1,'first');
start=max(onset-35,1);

lumps=findLumps(trials(N).t-trials(N).t(start),[trials(N).x trials(N).v],length(trials(N).t),0,4,1)

targ=trials(N).x(end,:);
for k=1:length(lumps)
    distrem(k)=norm(trials(N).x(lumps(k).inds(1),1:2)-targ);
    distcov(k)=norm(lumps(k).y(end,1:2)-lumps(k).y(1,1:2));
end

figure(8)
clf
plot(distrem,distcov,'-x')

%[distrem' ones(4,1)]\distcov'

%% Step 1: Guess at a peak and do serial decomposition

figure(1)
clf
hold on

t=trials(N).t;
lT=length(t);
y=[trials(N).x trials(N).v];
plot(1:lT,vecmag(y(:,3:4)));

%next=168; % 210 Pretty fractal breakdown
[blah,next]=max(vecmag(y(:,3:4)))
first=sum(vecmag(y(:,3:4)));
resid=y;

cycles=0;
while cycles<5
    cycles=cycles+1
    [lump,resid]=rulesFindLumps(t,resid,next);
    if length(lump.t)==1
        break
    end
    vm=vecmag(resid(:,3:4));
    c=rand(3,1);
    plot(1:lT,vm,'-','color',c);
    plot(lump.inds,vecmag(lump.y(:,3:4)),'.','color',c);
    [blah,next]=max(vm);
end
xlabel('Time, Samples (rate=200 Hz)')
ylabel('Speed, meters/sec')