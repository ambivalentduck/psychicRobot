clc
clear all

load(['./Data/christine_stiffness.mat']);

notouch=2360:3338;
held=3920:5170;

t=trials(1).t(notouch);
xnotouch=trials(1).x(notouch,1)-trials(1).x(notouch(1),1);
ynotouch=trials(1).x(notouch,2)-trials(1).x(notouch(1),2);

figure(1)
clf
plot(t,xnotouch,'b',t,ynotouch,'r')

std(xnotouch)