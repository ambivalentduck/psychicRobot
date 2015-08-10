clc
clear all
close all

load('free_exp_05stroke.mat')

ti=find((t>3.55)&(t<10));
srate=mean(gradient(t));

t=t(ti);
x=x(ti,:);
v=v(ti,:);
a=a(ti,:);
P=dot(a',v');
speed=vecmag(v);

figure(1)
clf
hold on
plot(t,speed,'b')

% A leading edge decomposition starts from an isolated peak and works
% forward one dot-product at a time.

f=find(t==3.7);

d=ones(length(t),1)*v(f,:);
plot(t,dot(v',d'),'r')
gd=.2*gradient(dot(v',d'))./gradient(t'); 
plot(t,gd,'k')
%To make this more efficient, calculate dot and delta dot in the while-loop

left=f;
while gd(left-1)>gd(left)
    left=left-1;
end
plot(t(left),gd(left),'gx')

right=f;
while gd(right+1)<gd(right)
    right=right+1;
end
plot(t(right),gd(right),'gx')
