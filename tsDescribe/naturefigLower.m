clc
clear all

load('../Data/curlkick/curlkick1Y.mat')

tcats=[trials.targetcat];
dcats=[trials.disturbcat];

figure(1)
clf
hold on

colors='rgb';
dirs=[1 2 4];

for DERP=1:3
f=find((tcats==dirs(DERP))&~dcats);
S=zeros(size(f));
RW=S;

for k=1:length(f)
    tr=trials(f(k));
    [S(k),RW(k)]=getTsMetric(tr.x,tr.v,tr.a,tr.t,tr.x(1,:),tr.x(end,:));
end

[f,x]=ecdf(S);
[shift,n,T]=fitShiftedGam(S)
plot(x,f,colors(DERP))
plot(x,gamcdf(x-shift,n,T),[colors(DERP),'--'])
end