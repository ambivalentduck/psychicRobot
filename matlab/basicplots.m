clc
clear all

load ../Data/1.mat
load ../Data/1extracted.mat

types=[trials.type];

STIFFNESS=2;

for k=1:3
    figure(k)
    clf
end

marker='-';
SPACE=2;
SATURATE=.7;

for k=1:length(trials)
   figure(trials(k).type)
   x=desiredTrajectories(k,STIFFNESS).xDesired(1:SPACE:end,1)-desiredTrajectories(k,STIFFNESS).xDesired(1,1);
   y=desiredTrajectories(k,STIFFNESS).xDesired(1:SPACE:end,2)-desiredTrajectories(k,STIFFNESS).xDesired(1,2);
   subplot(3,1,1)
   hold on
   if ~trials(k).long
       yoff=.1;
   else
       yoff=0;
   end
   plot(x,y+yoff,marker)
   pos=[trials(k).pos(1:SPACE:end,1)-trials(k).pos(1,1),trials(k).pos(1:SPACE:end,2)-trials(k).pos(1,2)+yoff];
   quiver(pos(:,1),pos(:,2),trials(k).force(1:SPACE:end,1),trials(k).force(1:SPACE:end,2),marker,'Color',[1, SATURATE SATURATE])
   subplot(3,1,2)
   hold on
   plot(x,vecmag(desiredTrajectories(k,STIFFNESS).vDesired(1:SPACE:end,:)),marker)
   plot(trials(k).pos(1:SPACE:end,1)-trials(k).pos(1,1),vecmag(trials(k).vel(1:SPACE:end,:)),marker,'Color',[1, SATURATE SATURATE])
   subplot(3,1,3)
   hold on
   plot(desiredTrajectories(k,STIFFNESS).time(1:SPACE:end),vecmag(desiredTrajectories(k,STIFFNESS).vDesired(1:SPACE:end,:)),marker)
   plot(trials(k).time(trials(k).first:SPACE:end)-trials(k).time(trials(k).first),vecmag(trials(k).vel(trials(k).first:SPACE:end,:)),marker,'Color',[1, SATURATE SATURATE])
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