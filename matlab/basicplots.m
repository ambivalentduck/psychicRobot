clc
clear all

load ../Data/1.mat
load ../Data/1extracted.mat

types=[trials.type];

STIFFNESS=1;

for k=1:3
    figure(k)
    clf
end

for k=1:length(trials)
   figure(trials(k).type)
   x=desiredTrajectories(k,STIFFNESS).xDesired(:,1)-desiredTrajectories(k,STIFFNESS).xDesired(1,1);
   y=desiredTrajectories(k,STIFFNESS).xDesired(:,2)-desiredTrajectories(k,STIFFNESS).xDesired(1,2);
   subplot(3,1,1)
   hold on
   plot(x,y)
   subplot(3,1,2)
   hold on
   plot(x,vecmag(desiredTrajectories(k,STIFFNESS).vDesired))
   subplot(3,1,3)
   hold on
   plot(desiredTrajectories(k,STIFFNESS).time,vecmag(desiredTrajectories(k,STIFFNESS).vDesired))
end

titles={'Early Kick: 10% completion','Late Kick: 50% completion','BLWN'};

for k=1:3
    figure(k)
    subplot(3,1,1)
    title(titles{k})
    axis equal
    ylabel('Robot y, meters')
    xlabel('Robot x, meters')
    subplot(3,1,2)
    ylabel('Speed, m/s')
    xlabel('Translated Robot x, meters')
    subplot(3,1,3)
    ylabel('Speed, m/s')
    xlabel('Time since onset of movement, seconds')
end