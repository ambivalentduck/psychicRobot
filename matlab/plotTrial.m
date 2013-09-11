function plotTrial(trial,params)

global kpgain

%f=find((vecmag(trial.accel)<.15)&(vecmag(trial.vel)<.001));
f=find(((trial.time-trial.time(1))<.5)'&(vecmag(trial.vel)<.005));

forcecalib=trial.force(f,:);
mforce=mean(forcecalib);
mq=mean(trial.qRobot(f,:)); %Hack that should be fair. We're not exactly moving.

force=trial.force;
mf=mforce';
for k=1:length(trial.time)
    q=sum(trial.qRobot(k,:)-mq);
    s=sin(q);
    c=cos(q);
    force(k,:)=force(k,:)-([c -s;s c]*mf)';
end
sum(sum(abs(force)))

%force=[trial.force(:,1)-mforce(1) trial.force(:,2)-mforce(2)];
t=trial.time-trial.time(1);

figure(1)
clf
subplot(3,1,1)
v=vecmag(trial.force);
N=6;
%fv=filter(ones(N,1)/N,1,v);
[vals,locs]=findpeaks(v,'minpeakheight',5);
locs=locs-N/2;
hold on
plot(trial.time,v,'b')
plot(trial.time(locs),v(locs-3),'rx')
w=f(end):locs(1)-2;
plot(trial.time(w),v(w),'k')

subplot(3,1,2)
hold on
vmf=vecmag(force);
plot(t,vmf)
w=find((vmf<3)&(vmf>.2)); %areas where force is middling

subplot(3,1,3)
X=[trial.qddotRobot trial.qdotRobot cos(trial.qdotRobot(:,2)).*(2*trial.qdotRobot(:,1)+trial.qdotRobot(:,2)) trial.qdotRobot(:,1).^2];

mtorque=mean(trial.torqueRobot(f,:));

for k=1:length(t)
    torque(k,:)=(robotfJ(trial.qRobot(k,:))*force(k,:)')';
end

y=torque;
Xw=X(w,:);
m=Xw\y(w,:)
plot(t,y,t,X*m,'.')

for k=1:length(t) %Comment line below to skip accounting for robot dynamics
    force(k,:)=force(k,:)-(robotfJ(trial.qRobot(k,:))\(X(k,:)*m)')';
end

subplot(3,1,2)
plot(t,vecmag(force),'r.')

xvaf=[trial.pos trial.vel trial.accel force];


figure(trial.rawnum)
clf
subplot(2,1,1)
hold on
plot(trial.pos(:,1),trial.pos(:,2),'b-o','MarkerSize',2)
quiver(trial.pos(:,1),trial.pos(:,2),trial.force(:,1),trial.force(:,2),'Color',[.5 .5 .5])

subplot(2,2,3)
hold on
v=vecmag(trial.vel);
plot(trial.pos(:,1),v,'b')
plot(trial.pos(f,1),v(f),'b.')
subplot(2,2,4)
hold on
plot(t,v,'b')
plot(t(f),v(f),'b.')

[trash,i]=min(abs(t-.96));

kpgain=1;
y=extract(t,xvaf,params,@armdynamics_inverted);
subplot(2,1,1)
plot(y(:,1),y(:,2),'k-o','MarkerSize',2)
subplot(2,2,3)
plot(y(:,1),vecmag(y(:,[3 4])),'k')
subplot(2,2,4)
plot(t,vecmag(y(:,[3 4])),'k')
plot(t,vecmag(y(:,[3 4])-trial.vel),'k-.')


y=extract(t,xvaf,params,@armdynamicsInvertedBurdet);
subplot(2,1,1)
plot(y(:,1),y(:,2),'r-o','MarkerSize',2)
subplot(2,2,3)
plot(y(:,1),vecmag(y(:,[3 4])),'r')
subplot(2,2,4)
hold on
plot(t,vecmag(y(:,[3 4])),'r')

for kpgain=.35 %linspace(.15,5,20)
y=extract(t,xvaf,params,'reflex');
subplot(2,1,1)
plot(y(:,1),y(:,2),'m-o','MarkerSize',2)
subplot(2,2,3)
plot(y(:,1),vecmag(y(:,[3 4])),'m')
subplot(2,2,4)
hold on
plot(t,vecmag(y(:,[3 4])),'m')
end

% kpgain=2;
% y=extract(t(i:end),xvaf(i:end,:),params,@armdynamics_inverted,y(i,:));
% subplot(2,1,1)
% plot(y(:,1),y(:,2),'m--v','MarkerSize',2)
% subplot(2,2,3)
% plot(y(:,1),vecmag(y(:,[3 4])),'m--')
% subplot(2,2,4)
% plot(t(i:end),vecmag(y(:,[3 4])),'m--')

subplot(2,1,1)
axis equal
legend('Measured','Measured Forces','Const Impedance','Torque-Scaled Impedance','Scaled Imp + Reflexes')
title(['Trial Number ',num2str(trial.rawnum)])
plot([trial.origin(1) trial.target(1)],[trial.origin(2) trial.target(2)],'Color',[.5 .5 .5])

subplot(2,2,3)
ylabel('Speed, m/s')
xlabel('Robot x, m')

subplot(2,2,4)
ylabel('Speed, m/s')
xlabel('Time, s')