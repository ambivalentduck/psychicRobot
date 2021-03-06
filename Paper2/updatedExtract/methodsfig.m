clc
clear all

for k=1:5
    P(k,:)=[96*(k-1)+1,96*k];
end

figure(1)
clf
hold on

plot([P(1,1),P(5,2)],0*[1 1],'k','linewidth',5)
plot([P(2,1),P(4,2)],1*[1 1],'k','linewidth',5)
plot([P(3,1),P(3,2)],3*[1 1],'k','linewidth',5)
plot([P(1,1),P(2,2)],2*[1 1],'k','linewidth',5)
plot([P(4,1),P(5,2)],2*[1 1],'k','linewidth',5)

xlim([0 5*96])
ylim([-.1 3.1])
set(gca,'ytick',[0 1 2 3])
set(gca,'yticklabels',{'Targetted Reaching in Three Directions','Turbulent Forces Applied','Cursor = Hand','Cursor = Intent'})

matlabfrag('figures/fig3raw')