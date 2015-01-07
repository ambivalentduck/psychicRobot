clc
clear all

load('3dmetrics.mat')

h=load('healthy.mat');
h=h.outs;
s=load('stroke.mat');
s=s.outs;

figure(1)
clf
subplot(1,2,1)
hold on
plot(zeros(size(h,1),1),h(:,1),'r.')
plot(ones(size(s,1),1),s(:,1),'r.')
xlim([-.5 1.5])
set(gca,'xtick',[0 1])
set(gca,'xticklabel',{'Healthy','Stroke'})

[yn,p]=ttest2(h(:,1),s(:,1))
title(['R^2 Boltzmann Distribution, T-test P=',num2str(p)])

subplot(1,2,2)
hold on
plot(zeros(size(h,1),1),h(:,3),'b.')
plot(ones(size(s,1),1),s(:,3),'b.')
xlim([-.5 1.5])
set(gca,'xtick',[0 1])
set(gca,'xticklabel',{'Healthy','Stroke'})

[yn,p]=ttest2(h(:,3),s(:,3))
title(['Boltzmann Temperature, T-test P=',num2str(p)])

set(gcf,'position',[76 11 1195 925])
print('-dtiff','-r300','summaryAll.tiff')
print('-dpng','-r300','summaryAll.png')

figure(2)
clf
subplot(1,2,1)
hold on
plot(zeros(size(h,1),1),h(:,2),'r.')
plot(ones(size(s,1),1),s(:,2),'r.')
xlim([-.5 1.5])
set(gca,'xtick',[0 1])
set(gca,'xticklabel',{'Healthy','Stroke'})

[yn,p]=ttest2(h(:,2),s(:,2))
title(['R^2 Boltzmann Distribution for Power (dot(a,v)), T-test P=',num2str(p)])

subplot(1,2,2)
hold on
plot(zeros(size(h,1),1),h(:,4),'b.')
plot(ones(size(s,1),1),s(:,4),'b.')
xlim([-.5 1.5])
set(gca,'xtick',[0 1])
set(gca,'xticklabel',{'Healthy','Stroke'})

[yn,p]=ttest2(h(:,4),s(:,4))
title(['Characteristic Power (essentially Wattage), T-test P=',num2str(p)])

set(gcf,'position',[76 11 1195 925])
print('-dtiff','-r300','summaryAllpower.tiff')
print('-dpng','-r300','summaryAllpower.png')