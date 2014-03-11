clc
clear all

colors=[1 0 0; 0 1 0; 0 0 1; 0 0 0];
U=7

figure(6)
    clf
for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    subplot(2,1,1)
    hold on
    plot(-1:U,W(:,1),'-','color',colors(k,:))
    ylabel('Fit Feedback Gain, unitless')
    set(gca,'xticklabel',[])
    subplot(2,1,2)
    hold on
    plot(-1:U,W(:,2),'.','color',colors(k,:))
    plot(-1:U,W(:,3),'-o','color',colors(k,:))
    ylabel('Mean Error, Nm')
    xlabel('Trials Used to Determine Weight')
    set(gca,'xticklabel',labels)
end

subplot(2,1,1)
ylim([0 1])