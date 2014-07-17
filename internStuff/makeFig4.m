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

units={'Max Error, cm','Stiffness, N/m','Error Span, cm'};
spacer=[.2 450 .14];

colSpan=.5;
subWidth=colSpan/10;
SUBS=8;

hand(:,2)=force(:,2)./hand(:,2);
cursor(:,2)=force(:,2)./cursor(:,2);

summaryMet=@median;

% The figure that will display all three plots on one figure
figure(1)
clf
set(gcf,'color','w')

map=[1 3 2];

for N=1:3
    
    for S = 1:SUBS
        col=[.7*rand(1,3)];
        dcol=1-col;
        
        subplot(3,1,map(N));
        hold on
        
        xoff=colSpan*((S-1)/(SUBS-1))-colSpan/2;
        
        for P=1:5
            if N~=2
                y=cursor(phases(:,P),S,N);
            else
                y=hand(phases(:,P),S,N);
            end
            
            plot(ones(96,1)*P+xoff+subWidth*(rand(96,1)-.5),y,'.','markersize',.00001,'color',col)
        end
        if N == 2
            xlabel('Trials')
        end
        
        ylabel(units{N})
        
    end
end
