clc
clear all
load('pulse1.mat');
close all

N=14;
onset=find(vecmag(trials(N).v)>.05,1,'first');
start=max(onset-35,1);

% lumps=findLumps(trials(N).t-trials(N).t(start),[trials(N).x trials(N).v],length(trials(N).t),0,4,1)
% 
% targ=trials(N).x(end,:);
% for k=1:length(lumps)
%     distrem(k)=norm(trials(N).x(lumps(k).inds(1),1:2)-targ);
%     distcov(k)=norm(lumps(k).y(end,1:2)-lumps(k).y(1,1:2));
% end
% 
% figure(8)
% clf
% plot(distrem,distcov,'-x')

%[distrem' ones(4,1)]\distcov'

%% Step 1: Guess at a peak and do serial decomposition

figure(1)
clf
hold on

t=trials(N).t;
lT=length(t);
y=[trials(N).x trials(N).v];

subplot(2,1,1)
hold on
plot(trials(N).x(:,1), trials(N).x(:,2),'k')

subplot(2,1,2)
hold on
plot(1:lT,vecmag(y(:,3:4)));

%next=371; % 210 Pretty fractal breakdown
[blah,next]=max(vecmag(y(:,3:4)))
first=sum(vecmag(y(:,3:4)));
resid=y;

tic
cycles=0;
while cycles<10
    cycles=cycles+1
    [lump,resid]=rulesFindLumps(t,resid,next);
    if length(lump.t)==1
        break
    end
    vm=vecmag(resid(:,3:4));
    c=rand(3,1);
    subplot(2,1,1)
    plot(lump.y(:,1),lump.y(:,2),'.','color',c)
    subplot(2,1,2)
    plot(1:lT,vm,'-','color',c);
    plot(lump.inds,vecmag(lump.y(:,3:4)),'.','color',c);
    [blah,next]=max(vm);
end
toc
xlabel('Time, Samples (rate=200 Hz)')
ylabel('Speed, meters/sec')

%% Step 2: Try "interactive" leading edge decomposition

figure(2)
clf
hold on

t=trials(N).t;
lT=length(t);
y=trials(N).v;
resid=y;

inds=[143 213 273];

plot(1:lT,vecmag(y));

for k=1:length(inds)-1
    z=zeros(size(resid));
    z(inds(k):inds(k+1),:)=resid(inds(k):inds(k+1),:);
    z(inds(k+1)+1:2*inds(k+1)-inds(k),:)=resid(inds(k+1)-1:-1:inds(k),:);
    vm=vecmag(resid);
    c=rand(3,1);
    resid=resid-z;
    plot(1:lT,vecmag(resid),'-','color',c);
    plot(1:lT,vecmag(z),'.','color',c);
end

xlabel('Time, Samples (rate=200 Hz)')
ylabel('Speed, meters/sec')