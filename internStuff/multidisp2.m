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

hand(:,2)=force(:,2)./hand(:,2);
cursor(:,2)=force(:,2)./cursor(:,2);

summaryMet=@median;

for N=1:3
    figure(N)
    clf
    hold on


    for S = 1:8
        col=[rand(1,2) 0];
        dcol=[col(1:2) 1];

        %plot(allT,cursor(:,N)+S*spacer(N),'.','color',col)
        plot([2*96+1,3*96],summaryMet(cursor(phases(:,3),N))+[0 0]+S*spacer(N),'-','color',col)
        %plot(allT(phases(:,3)),hand(phases(:,3),N)+S*spacer(N),'x','color',dcol)
        plot([2*96+1,3*96],summaryMet(hand(phases(:,3),N))+[0 0]+S*spacer(N),'-','color',dcol)
        plot([97 4*96],summaryMet(cursor(phases(:,2)|phases(:,4),N))+[0 0]+S*spacer(N),'-.','color','r')

        %    errave = (mean(errors(phase2(curs0),1))+mean(errors(phase4(curs0),1)))/2;
        %    plot([97 4*96],errave*[1 1]+num*.2,'--','color','red')
        %axis ([0 5*96 0 .2])
    end
    title(titles{N})
    for g=1:4
        line([g*96,g*96],ylim)
    end
end


