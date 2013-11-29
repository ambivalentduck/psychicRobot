clc
clear all

t=(0:.01:1)';

y=[t cos(pi*t) [0; ones(length(t)-2,1);0] -pi*sin(pi*t)];

o10=ones(10,1);
y=[o10*y(1,:)
    y;
    o10*y(end,:)];
    
t=[.01*(-10:-1)'; t; 1+.01*(1:10)'];

figure(1)
clf
hold on
plot(y(:,1),y(:,2))

findLumps(t,y,length(t),0,10,2)