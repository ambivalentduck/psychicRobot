clc
clear all
load('pulse1.mat');

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