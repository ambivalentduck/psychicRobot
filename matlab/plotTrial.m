function plotTrial(trial,params)

figure(trial.rawnum)
clf
subplot(2,1,1)

forcecalib=trial.force((((trial.time-trial.time(1))<.3)&(vecmag(trial.vel)<.01)),:);
mforce=mean(forcecalib);
force=[trial.force(:,1)-mforce(1) trial.force(:,2)-mforce(2)];
xvaf=[trial.pos trial.vel trial.accel -force];
t=trial.time;
t=linspace(t(1),t(end),length(t));
y=extract(t,xvaf,params,@armdynamicsInvertedBurdet);
%y=extract(t,xvaf,params,@armdynamics_inverted);

hold on
plot(trial.pos(:,1),trial.pos(:,2))
plot(y(:,1),y(:,2),'r')
axis equal

subplot(2,1,2)
% distreal=cumsum(sqrt(gradient(trial.pos(:,1)).^2+gradient(trial.pos(:,2)).^2));
% distext=cumsum(sqrt(gradient(y(:,1)).^2+gradient(y(:,2)).^2));
plot(trial.pos(:,1),vecmag(trial.vel),'b',y(:,1),vecmag(y(:,[3 4])),'r')