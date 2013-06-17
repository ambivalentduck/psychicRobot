clc
clear all

load ../Data/3.mat
load ../Data/3extracted.mat

types=[trials.type];

STIFFNESS=3;

for k=1:3
    figure(k)
    clf
end

marker='-';
SPACE=2;
SATURATE=.7;

for k=1:length(trials)
   figure(trials(k).type)
   xoff=-trials(k).pos(1,1);
   yoff=trials(k).updown/50;
   if ~trials(k).long
       xoff=xoff+1;
   end
   x=desiredTrajectories(k,STIFFNESS).xDesired(1:SPACE:end,1)-xoff;
   y=desiredTrajectories(k,STIFFNESS).xDesired(1:SPACE:end,2)-desiredTrajectories(k,STIFFNESS).xDesired(1,2);
   subplot(5,1,1:3)
   hold on
   plot(x,y+yoff,marker)
   plot(x(1),y(1)+yoff,'rx')
   pos=[trials(k).pos(1:SPACE:end,1)-xoff,trials(k).pos(1:SPACE:end,2)-trials(k).pos(1,2)+yoff];
   plot(pos(:,1),pos(:,2),marker,'Color',[1, SATURATE SATURATE])
   quiver(pos(:,1),pos(:,2),trials(k).force(1:SPACE:end,1),trials(k).force(1:SPACE:end,2),marker,'Color',[1, SATURATE SATURATE])
   subplot(5,1,4)
   hold on
   plot(x,vecmag(desiredTrajectories(k,STIFFNESS).vDesired(1:SPACE:end,:)),marker)
   plot(trials(k).pos(1:SPACE:end,1)-xoff,vecmag(trials(k).vel(1:SPACE:end,:)),marker,'Color',[1, SATURATE SATURATE])
   subplot(5,1,5)
   hold on
   plot(desiredTrajectories(k,STIFFNESS).time(1:SPACE:end),vecmag(desiredTrajectories(k,STIFFNESS).vDesired(1:SPACE:end,:)),marker)
   plot(trials(k).time(trials(k).first:SPACE:end)-trials(k).time(trials(k).first),vecmag(trials(k).vel(trials(k).first:SPACE:end,:)),marker,'Color',[1, SATURATE SATURATE])
end

titles={'Early Kick: 10% completion','Late Kick: 50% completion','BLWN'};

for k=1:3
    figure(k)
    subplot(5,1,1:3)
    title(titles{k})
    axis equal
    ylabel('Robot y, meters')
    xlabel('Robot x, meters')
    subplot(5,1,4)
    ylabel('Speed, m/s')
    xlabel('Translated Robot x, meters')
    subplot(5,1,5)
    ylabel('Speed, m/s')
    xlabel('Time since onset of movement, seconds')
end