clc
clear all

load ../Data/11.mat
load ../Data/11extracted.mat

types=[trials.type];

STIFFNESS=1;

for k=1:3
    figure(k)
    clf
end

marker='-';
SPACE=2;
SATURATE=.7;

for k=1:length(trials)
   figure(trials(k).type)
   xoff=trials(k).pos(1,1)+.4*sign(trials(k).updown);
   yoff=0; %trials(k).updown/50;
   if ~trials(k).long
       yoff=yoff+.25;
   end
   %subplot(5,1,1:3)
   hold on
   
   xshad=trials(k).des(1:SPACE:end,1);
   yshad=trials(k).des(1:SPACE:end,2)-trials(k).des(1,2);
   plot(xshad-xoff,yshad+yoff,marker,'Color',[SATURATE 1 SATURATE])
   plot(xshad(1)-xoff,yshad(1)+yoff,'rx')

   xburdet=desiredTrajectories(k,1).xDesired(:,1);
   yburdet=desiredTrajectories(k,1).xDesired(:,2)-desiredTrajectories(k,1).xDesired(1,2);
   plot(xburdet-xoff,yburdet+yoff,'.','Color',[SATURATE SATURATE 1])
   
   xburdetreflexes=desiredTrajectories(k,2).xDesired(:,1);
   yburdetreflexes=desiredTrajectories(k,2).xDesired(:,2)-desiredTrajectories(k,2).xDesired(1,2);
   plot(xburdet-xoff,yburdet+yoff,marker,'Color',[1-SATURATE 1-SATURATE 1-SATURATE])
   
   pos=[trials(k).pos(1:SPACE:end,1)-xoff,trials(k).pos(1:SPACE:end,2)-trials(k).pos(1,2)+yoff];
   plot(pos(:,1),pos(:,2),marker,'Color',[1, SATURATE SATURATE])
   quiver(pos(:,1),pos(:,2),trials(k).force(1:SPACE:end,1),trials(k).force(1:SPACE:end,2),marker,'Color',[1, SATURATE SATURATE])
   plot([trials(k).origin(1) trials(k).target(1)]-xoff,[0 0]+yoff,'m-')
   
%    subplot(5,1,4)
%    hold on
%    plot(x,vecmag(desiredTrajectories(k,STIFFNESS).vDesired(1:SPACE:end,:)),marker)
%    plot(trials(k).pos(1:SPACE:end,1)-xoff,vecmag(trials(k).vel(1:SPACE:end,:)),marker,'Color',[1, SATURATE SATURATE])
%    subplot(5,1,5)
%    hold on
%    plot(desiredTrajectories(k,STIFFNESS).time(1:SPACE:end),vecmag(desiredTrajectories(k,STIFFNESS).vDesired(1:SPACE:end,:)),marker)
%    plot(trials(k).time(trials(k).first:SPACE:end)-trials(k).time(trials(k).first),vecmag(trials(k).vel(trials(k).first:SPACE:end,:)),marker,'Color',[1, SATURATE SATURATE])
end

titles={'Early Kick: 10% completion','Late Kick: 50% completion','BLWN'};

% for k=1:3
%     figure(k)
%     subplot(5,1,1:3)
%     title(titles{k})
%     axis equal
%     ylabel('Robot y, meters')
%     xlabel('Robot x, meters')
%     subplot(5,1,4)
%     ylabel('Speed, m/s')
%     xlabel('Translated Robot x, meters')
%     subplot(5,1,5)
%     ylabel('Speed, m/s')
%     xlabel('Time since onset of movement, seconds')
% end