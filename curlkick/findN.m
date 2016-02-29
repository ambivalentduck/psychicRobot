clc
clear all

load ../Data/curlkick/curlkick1Y.mat

doPlots=0;
clean_init=1;

f=find(([trials.targetcat]~=0)&(~[trials.disturbcat]));

Ts=zeros(size(f));
R=zeros(size(f));

for ff=1:length(f)
    T=f(ff);
    [Ts(ff),R]=getRmetrics(trials(ff).x,trials(ff).v,trials(ff).a,trials(ff).t,trials(ff).x(1,:),trials(ff).x(end,:));
end

figure(3)
clf
[shift,n,T,rmse]=fitShiftedGam(Ts.^-2,1)
