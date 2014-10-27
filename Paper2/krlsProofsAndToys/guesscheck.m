clc
clear all
close all

t=0:.01:.5;

[kern,fop]=slmj5op(t,.25,.5)

fop=fop*.1;

figure(1)
clf
hold on
plot(t,fop)
plot(t,fop(:,2).*fop(:,3))

figure(2)
clf
hist(abs(fop(:,2).*fop(:,3)),40)

lower=tc-ts/2;
upper=tc+ts/2;

% ta=(t-lower)/ts;
% out=[10*ta.^3-15*ta.^4+6*ta.^5,
%     (30*ta.^2-60*ta.^3+30*ta.^4)/ts,
%     (60*ta-180*ta.^2+120*ta.^3)/(ts^2)];