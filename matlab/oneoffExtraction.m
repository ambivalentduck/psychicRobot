clc
clear all
load ../Data/729.mat
close all

global kpgain
kpgain=.25;

for k=1:50
    %continue
    if trials(k).rawnum<11
        figure(5001)
        if trials(k).pos(1,1)<0
            subplot(2,2,2)
        else
            subplot(2,2,1)
        end
        hold on
        plot(trials(k).pos(:,1),trials(k).pos(:,2),'g')
        axis equal
        figure(5001)
        if trials(k).pos(1,1)<0
            subplot(2,2,4)
        else
            subplot(2,2,3)
        end
        hold on
        plot(trials(k).pos(:,1),trials(k).pos(:,2),'g')
        axis equal
        
        
        figure(5002)
        if trials(k).pos(1,1)<0
            subplot(2,1,1)
        else
            subplot(2,1,2)
        end
        hold on
        plot(trials(k).pos(:,1),trials(k).pos(:,2),'g')
        axis equal
        continue
    end
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
    xvaf=[trials(k).pos trials(k).vel trials(k).accel force];
    %y=extract(t,xvaf,params,'reflex');
    y=extract(t,xvaf,params,@armdynamicsInvertedBurdet);
    %y=extract(t,xvaf,params,@armdynamics_inverted);
    %     if trials(k).rawnum>=30
    %         figure(trials(k).rawnum)
    %         clf
    %         subplot(2,1,1)
    %         hold on
    %         plot(trials(k).pos(:,1),trials(k).pos(:,2))
    %         plot(y(:,1),y(:,2),'r')
    %         axis equal
    %         subplot(2,1,2)
    %         plot(t,vecmag(y(:,3:4)))
    %     end


    if trials(k).pos(1,1)<0
        left=1;
    else
        left=0;
    end

    if trials(k).white==0
        figure(5001)
        early=trials(k).early~=0;
        subplot(2,2,early*2+left+1)
    else
        figure(5002)
        subplot(2,1,left+1)
    end

    hold on
    plot(trials(k).pos(:,1),trials(k).pos(:,2),'b-o','MarkerSize',2)
    plot(y(:,1),y(:,2),'r-o','MarkerSize',2)
    try
        plot(trials(k).pos(still(end)+50,1),trials(k).pos(still(end)+50,2),'bx')
        plot(y(still(end)+50,1),y(still(end)+50,2),'rx')
    catch
        k
    end
    axis equal
end

shapes=[trials(51:66).shape];
us=unique(shapes);
lus=length(us);
mags=[trials(51:66).white];
u=unique(mags);
lu1=length(u)+1;

% figure(20)
% clf
% for k=51:66
% 
%     mag=trials(k).white;
%     f=find(u==mag);
%     if trials(k).vision==0
%         f=f+1;
%     end
%     subplot(lus,lu1,trials(k).shape+1+(f-1)*lus)
%     hold on
% 
%     still=find(((trials(k).time-trials(k).time(1))<.3)&(vecmag(trials(k).vel)<.01));
%     forcecalib=trials(k).force(still,:);
%     mforce=mean(forcecalib);
%     force=[trials(k).force(:,1)-mforce(1) trials(k).force(:,2)-mforce(2)];
%     xvaf=[trials(k).pos trials(k).vel trials(k).accel -force];
%     t=trials(k).time;
%     t=linspace(t(1),t(end),length(t));
%     %y=extract(t,xvaf,params,@armdynamicsInvertedBurdet);
%     y=extract(t,xvaf,params,@armdynamics_inverted);
% 
%     plot(trials(k).pos(:,1),trials(k).pos(:,2))
%     plot(y(:,1),y(:,2),'r')
%     axis equal
%     if trials(k).shape==0
%         switch f
%             case 1
%                 ylabel('Vision')
%             case 2
%                 ylabel('No Vision')
%             case 3
%                 ylabel('No Vision + 1 BLWN')
%             case 4
%                 ylabel('No Vision + 1.5 BLWN')
%         end
%     end
% end

cleanup