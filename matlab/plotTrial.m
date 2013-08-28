function plotTrial(trial,params)

global kpgain

kpgain=.25;

figure(trial.rawnum)
clf
subplot(2,1,1)

forcecalib=trial.force((((trial.time-trial.time(1))<.3)&(vecmag(trial.vel)<.01)),:);
mforce=mean(forcecalib);
force=[trial.force(:,1)-mforce(1) trial.force(:,2)-mforce(2)];
xvaf=[trial.pos trial.vel trial.accel force];
t=trial.time;
t=linspace(t(1),t(end),length(t));

hold on
plot(trial.pos(:,1),trial.pos(:,2),'b-o','MarkerSize',2)
quiver(trial.pos(:,1),trial.pos(:,2),trial.force(:,1),trial.force(:,2),'Color',[.5 .5 .5])

subplot(2,2,3)
hold on
plot(trial.pos(:,1),vecmag(trial.vel),'b')
subplot(2,2,4)
hold on
plot(trial.time,vecmag(trial.vel),'b')


y=extract(t,xvaf,params,@armdynamics_inverted);
subplot(2,1,1)
plot(y(:,1),y(:,2),'k-o','MarkerSize',2)
subplot(2,2,3)
plot(y(:,1),vecmag(y(:,[3 4])),'k')
subplot(2,2,4)
hold on
plot(trial.time,vecmag(y(:,[3 4])),'k')


y=extract(t,xvaf,params,@armdynamicsInvertedBurdet);
subplot(2,1,1)
plot(y(:,1),y(:,2),'r-o','MarkerSize',2)
subplot(2,2,3)
plot(y(:,1),vecmag(y(:,[3 4])),'r')
subplot(2,2,4)
hold on
plot(trial.time,vecmag(y(:,[3 4])),'r')


y=extract(t,xvaf,params,'reflex');
subplot(2,1,1)
plot(y(:,1),y(:,2),'m-o','MarkerSize',2)
subplot(2,2,3)
plot(y(:,1),vecmag(y(:,[3 4])),'m')
subplot(2,2,4)
hold on
plot(trial.time,vecmag(y(:,[3 4])),'m')

subplot(2,1,1)
axis equal
legend('Measured','Measured Forces','Const Impedance','Torque-Scaled Impedance','Scaled Imp + Reflexes')
title(['Trial Number ',num2str(trial.rawnum)])

subplot(2,2,3)
ylabel('Speed, m/s')
xlabel('Robot x, m')

subplot(2,2,4)
ylabel('Speed, m/s')
xlabel('Time, s')