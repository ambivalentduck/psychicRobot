clc
clear all

S=1;
L=[.15 0];
t=0:.01:1;
x=-(cos(pi*t/S)-1)/2;
v=pi*sin(pi*t/S)/(2*S);

figure(17)
clf
hold on
plot(t,v,t,x)


[lumps,resid]=findLumps(t',v'*L,(1:length(t))')

