clc
clear all

load ../Data/10.mat
load ../Data/10extracted.mat

types=[trials.type];

STIFFNESS=1;

for k=1:4
    figure(k)
    clf
end

G=unique([trials.updown]);

marker='-';
SPACE=2;
SATURATE=.3;

for k=1:length(trials)
   figure(trials(k).shape+1)
   F=find(G==trials(k).updown);
   subplot(length(G),3,3*F)
   xoff=0;
   yoff=0;
   hold on
   pos=[trials(k).pos(1:SPACE:end,1)-xoff,trials(k).pos(1:SPACE:end,2)-trials(k).pos(1,2)+yoff];
   plot(pos(:,1),pos(:,2),'b-')
   %quiver(pos(:,1),pos(:,2),trials(k).force(1:SPACE:end,1),trials(k).force(1:SPACE:end,2),'b-')
   axis equal
   ylabel(num2str(G(F)))
end

load ../Data/9.mat
load ../Data/9extracted.mat

types=[trials.type];

STIFFNESS=1;

G=unique([trials.updown]);

marker='-';
SPACE=2;
SATURATE=.3;

for k=1:length(trials)
   figure(trials(k).shape+1)
   F=find(G==trials(k).updown);
   xoff=0;
   yoff=0;
   subplot(length(G),3,3*F-1)
   x=desiredTrajectories(k,STIFFNESS).xDesired(1:SPACE:end,1)-xoff;
   y=desiredTrajectories(k,STIFFNESS).xDesired(1:SPACE:end,2)-desiredTrajectories(k,STIFFNESS).xDesired(1,2);
   hold on
   plot(x,y+yoff,marker)
   pos=[trials(k).pos(1:SPACE:end,1)-xoff,trials(k).pos(1:SPACE:end,2)-trials(k).pos(1,2)+yoff];
      subplot(length(G),3,3*F-2)
      hold on
   plot(pos(:,1),pos(:,2),marker,'Color',[1, SATURATE SATURATE])
   %quiver(pos(:,1),pos(:,2),trials(k).force(1:SPACE:end,1),trials(k).force(1:SPACE:end,2),marker,'Color',[1, SATURATE SATURATE])
   axis equal
   ylabel(num2str(G(F)))
end

for k=1:4
    figure(k)
    subplot(5,3,1)
    title('Hand - No Vision')
    subplot(5,3,2)
    title('Extracted - No Vision')
    subplot(5,3,3)
    title('Extracted - Shown')
end

    
