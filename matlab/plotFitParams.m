clc
clear all

names={'lc1','lc2','m1','m2','I1','I2','Kp0_11','Kp0_12','Kp0_21','Kp0_22','Kp1_11','Kp1_12','Kp1_21','Kp1_22'};
kp0sm=[15, 6;6, 16];

fitp=zeros(8,14);
nomp=zeros(8,14);

if ~exist('params.mat','file')
    for S=1:8
        findWhiteKpGain(S)
        [fitp(S,:),nomp(S,:)]=fitAllParams(S);
    end
    save('params.mat','fitp','nomp')
end

figure(1)
clf
for k=1:14
    subplot(1,14,k)
    hold on
    plot(1:8,(fitp(:,k)-nomp(:,k))./nomp(:,k),'rx')
    plot([0 1],[0 0],'m-')
    yl=ylim;
    myl=max(abs(yl));
    ylim([-myl myl])
    title(names{k})
    if k==1
        ylabel('Fit-Nominal Discrepancy, normalized by nominal')
    end
    set(gca,'xtick',[])
    xlim([1 8])
end