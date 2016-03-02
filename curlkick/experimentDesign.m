clc
clear all

load ../Data/curlkick/curlkick1Y.mat

f=find(([trials.targetcat]~=0)&([trials.disturbcat]))
length(f)/length(trials)

figure(1)
clf
subplot(2,1,1)
plot(f,'.')
subplot(2,1,2)
ecdf(diff(f))

figure(2)
clf
hold on
for k=1:length(f)
    t=trials(f(k)).t;
    plot(t-t(1),vecmag(trials(f(k)).f))
end