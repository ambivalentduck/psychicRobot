clc
clear all
close all

load('../Data/11.mat'); %just for params
global kp kd

xs=[trials.target];
xs=xs(1:2:end-1);
u=unique(xs);
y=.5;

coeff=calcminjerk([u(1) y],[u(2) y],[0 0],[0 0],[0 0],[0 0],0,.5);
t=0:.005:.5;
[x,v,a]=minjerk(coeff,t);
x=x';
v=v';
a=a';

figure(1)
hold on
plot(x(:,1),x(:,2),'b.',[u(1) u(2)],[y y],'m')
axis equal

figure(2)
hold on
plot(t,vecmag(v),'b.')
axis equal

xvaf=[x v a 0*x];
kd=[2.3 .09; .09 2.4];
kp=SnMKpGainStruct(1.5);
kp=kp{1};
y=extract(t,xvaf,params);
figure(1)
plot(y(:,1),y(:,2),'ro')
figure(2)
plot(t,vecmag(y(:,3:4)),'ro')

xvaf(:,8)=xvaf(:,8)+5;
y=extract(t,xvaf,params);
figure(1)
plot(y(:,1),y(:,2),'g.')
figure(2)
plot(t,vecmag(y(:,3:4)),'g.')

x0=fkin(ikin(x(1,:))-kp\[0; 1]);
y=extract(t,xvaf,params,[x0' 0 0]);
figure(1)
plot(y(:,1),y(:,2),'k.')
figure(2)
plot(t,vecmag(y(:,3:4)),'k.')

xvaf=[x v a 0*x];
xvaf(:,7:8)=xvaf(:,7:8)+10*(rand(size(xvaf(:,7:8)))-.5);
y=extract(t,xvaf,params);
figure(1)
plot(y(:,1),y(:,2),'c.')
figure(2)
plot(t,vecmag(y(:,3:4)),'c.')


figure(1)
title('Surrogate Data')
legend('Intent','Rhumb Line','Extraction','Bias in Force Sensor, 1N','Bias, but starting at stead-state error','Uniform Random Force Error, span 10N')

figure(2)
title('Surrogate Data')
legend('Intent','Extraction','Bias in Force Sensor, 1N','Bias, but starting at stead-state error','Uniform Random Force Error, span 10N')
ylabel('Speed,m/s')
xlabel('Time, s')