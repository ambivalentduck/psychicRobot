clear;
clc;
load ./Data/output27.mat

T=105;
f=trials(T).f;
t=trials(T).t;

figure(1)
hold on
plot(t,f(:,1),'r')
plot(t,f(:,2))
xlabel('Time (s)')
ylabel('Force (N)')
set(gcf,'color','w')

legend('Force in x position','Force in y position')
