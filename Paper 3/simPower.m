clc
clear all

M=.85;

%Generate smooth forces

L=10000;
filtn=1*64;
filtType='loess';
tstep=.005;
t=0:tstep:tstep*L;
t=t';
f=.5*randn(length(t),2);
f=[smooth(f(:,1),filtn,filtType) smooth(f(:,2),filtn,filtType)];

figure(1)
clf
subplot(2,4,1)
hold on
plot(t,f(:,1),'b')
plot(t,f(:,2),'r')
xlabel('Time')
ylabel('Force')

a=f/M;
v=[cumsum(a(:,1)) cumsum(a(:,2))]*tstep;

subplot(2,4,2)
hold on
plot(t,v(:,1),'b')
plot(t,v(:,2),'r')
xlabel('Time')
ylabel('Velocity')

p=[cumsum(v(:,1)) cumsum(v(:,2))]*tstep;

subplot(2,4,3)
hold on
plot(p(:,1),p(:,2),'b')
ylabel('Position, y')
xlabel('Position, x')
axis equal

P=dot(v',M*a');
U=cumsum(P);

subplot(2,4,4)
hold on
plot(t,P,'b')
xlabel('Time')
ylabel('Power')

subplot(2,4,5)
hold on
plot(t,U,'b')
xlabel('Time')
ylabel('Potential')

subplot(2,4,6)
hold on
nbins=50;
[counts,bins]=hist(abs(P),nbins);
bar(bins,counts)
nzbins=bins(counts>0);
nzcounts=counts(counts>0);
logcounts=log(nzcounts);
Wz=[logcounts' ones(length(nzcounts),1)]\nzbins';
plot(bins,exp((bins-Wz(2))/Wz(1)),'r')
R=cov(nzbins,logcounts)/(std(nzbins)*std(logcounts));
R2=R(1,2)^2;
title(['R^2=',num2str(R2)])