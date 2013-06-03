clc
clear all
close all

out=genNewExp;
trials=length(out)
time=3*trials/60

figure(1)
x=[-.21 .26];
y=[.37 .6];
hold on
plot([x(1) x(1) x(2) x(2) x(1)],[y(1) y(2) y(2) y(1) y(1)],'r')
plot(out(:,1),out(:,2),'b-x')
axis equal

figure(2)
s=sum(out(:,3:5)~=0,2)>0;
f=find(s&(out(:,6)==-1));
catches=length(f)
f=f(2:end)-f(1:end-1);
hist(f,1:24)
    