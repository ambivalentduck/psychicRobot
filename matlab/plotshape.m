clc
clear all
close all

load ../Data/12.mat
load ../Data/12extracted.mat

types=[trials.type];

for k=1:4
    figure(k)
    clf
end

G=unique([trials(1:88).white]);

marker='-';
SPACE=2;
SATURATE=.3;

for k=1:87
   %figure(k)
   figure(trials(k+1).shape+1)
   F=find(G==trials(k+1).updown);
   xoff=0;
   yoff=0;

   subplot(length(G),4,4*F-3)
   hold on
   plot(trials(k).pos(1:SPACE:end,1),trials(k).pos(1:SPACE:end,2),'b-')
   %quiver(pos(:,1),pos(:,2),trials(k).force(1:SPACE:end,1),trials(k).force(1:SPACE:end,2),'b-')
   axis equal
   ylabel(num2str(G(F)))
   
   subplot(length(G),4,4*F-2)
   hold on
   plot(trials(k).des(1:SPACE:end,1),trials(k).des(1:SPACE:end,2),'b-')
   axis equal
   
   subplot(length(G),4,4*F-1)
   hold on
   plot(desiredTrajectories(k,1).xDesired(:,1),desiredTrajectories(k,1).xDesired(:,2),'b-')
   axis equal
   
   subplot(length(G),4,4*F)
   hold on
   plot(desiredTrajectories(k,2).xDesired(:,1),desiredTrajectories(k,2).xDesired(:,2),'b-')
   axis equal
end


for k=1:4
    figure(k)
    subplot(5,4,1)
    title('Hand')
    subplot(5,4,2)
    title('Extracted Shad')
    subplot(5,4,3)
    title('Extracted Burdet')
    subplot(5,4,4)
    title('Extracted BurdetReflexes')
end

    
