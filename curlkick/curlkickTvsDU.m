clc
clear all

load ../Data/curlkick/curlkick1g.mat

f=find(([trials.targetcat]==1)&(~[trials.disturbcat]))

trials(f)
ts=0*f;
rw=0*f;

for k=1:length(f)
    tfk=trials(f(k));
    [ts(k),rw(k)]=getTsMetric(tfk.x,tfk.v,tfk.a,tfk.t,tfk.x(1,:),tfk.x(end,:));
end

ts_orig=ts;
%ts=ts.^-2

[f,x]=ecdf(ts);
figure(1)
clf
hold on
plot(x,f)
p=gamfit(ts)
plot(x,gamcdf(x,p(1),p(2)),'r')