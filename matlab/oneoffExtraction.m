clc
clear all
load ../Data/500.mat
close all

global kp kd
kd=[2.3 .09; .09 2.4];
kp=SnMKpGainStruct(1);
kp=kp{1};
kp=[10.8 2.83; 2.51 8.67];

for k=1:length(trials)-1
    t=trials(k).time;
    t=linspace(t(1),t(end),length(t));
    stillL=find((vecmag([trials(k).pos(:,1)-.15,trials(k).pos(:,2)-.5])<.005)&(vecmag(trials(k).vel)<.01));
    stillR=find((vecmag([trials(k).pos(:,1)+.15,trials(k).pos(:,2)-.5])<.005)&(vecmag(trials(k).vel)<.01));
    if length(stillL)<length(stillR)
        still=stillR;
    else
        still=stillL;
    end
    forcecalib=trials(k).force(still,:);
    mforce=mean(forcecalib);
    force=[trials(k).force(:,1)-mforce(1) trials(k).force(:,2)-mforce(2)];
    xvaf=[trials(k).pos trials(k).vel trials(k).accel -force];
    y=extract(t,xvaf,params);
%     figure(trials(k).rawnum)
%     clf
%     subplot(2,1,1)
%     plot(y(:,1),y(:,2))
%     axis equal
%     subplot(2,1,2)
%     plot(t,vecmag(y(:,3:4)))

    if trials(k).rawnum<11
        figure(5000)
    elseif trials(k).rawnum<31
        figure(5001)
    else
        figure(5002)
    end
    
    if trials(k).pos(1,1)<0
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    
    hold on
    plot(trials(k).pos(:,1),trials(k).pos(:,2))
    plot(y(:,1),y(:,2),'r')
    axis equal
end