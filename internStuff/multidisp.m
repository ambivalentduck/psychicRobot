% This is a function that takes in an output file and plots error(on the y axis) against
% trial(on the x axis) for baseline, extraction, and the secondary
% baseline(washout) all next to each other for all 8 patients
clc
clear all
close all

load('./Data/errors.mat')

allT=1:5*96;
phases=zeros(5*96,5);
for k=1:5
    phases(96*(k-1)+1:k*96,k)=1;
end
phases=phases==1;

% builds one display for all of the folloing data

titles={'Max Perp','RMS','Straightness'};
spacer=[.2 450 .14];

SHIFT=.3;

hand(:,:,2)=force(:,:,2)./hand(:,:,2);
cursor(:,:,2)=force(:,:,2)./cursor(:,:,2);

summaryMet=@median;

% The figure that will display all three plots on one figure 
figure(1) 
% sets the background to white
set(gcf,'color','w')
clf
for N=1:3
%     figure(N)
%     clf
%     hold on


    for S = 1:8
        col=[0 .5 .5]; %[rand(1,2) 0];
        dcol=col %[col(1:2) 1];
        
%         subplot iteration step
        subplot(1,3,N);
        
%         Plot what you need
        hold on
        plot(allT,cursor(:,S,N)+S*spacer(N),'.','color',col,'MarkerSize',.01)
        plot([2*96+1,3*96],summaryMet(cursor(phases(:,3),S,N))+[0 0]+S*spacer(N),'-','color','m','linewidth',2)
        plot(allT(phases(:,3)),hand(phases(:,3),S,N)+S*spacer(N),'x','color',dcol)
        plot([2*96+1,3*96],summaryMet(hand(phases(:,3),S,N))+[0 0]+S*spacer(N),'-','color','k', 'linewidth',2)
        plot([97 4*96],summaryMet(cursor(phases(:,2)|phases(:,4),S,N))+[0 0]+S*spacer(N),'-.','color','r','linewidth',2)
        
%         sets the appropriate labels depending on the current graph
        if N == 1
%             Max perpendicular distance
            xlabel('Trials')
            ylabel('Meters')
        elseif N == 2
%             RMS
            xlabel('Trials')
            ylabel('Newtons per Meter')
        elseif N == 3
%             Straightness 
            xlabel('Trials')
            ylabel('Meters')
        end 
        
        
        
        %    errave = (mean(errors(phase2(curs0),1))+mean(errors(phase4(curs0),1)))/2;
        %    plot([97 4*96],errave*[1 1]+num*.2,'--','color','red')
        %axis ([0 5*96 0 .2])
    end
    title(titles{N})
    
%     Graphs the verticle lines that separate the 5 phases 
    for g=1:4
        line([g*96,g*96],ylim)
    end
end


