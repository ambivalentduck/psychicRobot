clc
clear all

load stroke.mat
stroke=outs;

load healthy.mat
healthy=outs;

figure(1)
clf

fields={'W','n','T'};

nsp=length(fields);

for k=1:length(fields)
    subplot(1,nsp,k)
    nh=[healthy.(fields{k})];
    ns=[stroke.(fields{k})];
    plot([ones(length(nh),1); 2*ones(length(ns),1)],[nh ns]','.')
    xlim([.5 2.5])
    set(gca,'xtick',[1 2])
    set(gca,'xticklabels',{'Healthy','Stroke'})
    title(regexprep(fields(k),'_',' '))
end

figure(2)
clf

fields={'power_symmetry','power_rmse','speed_rmse'};

nsp=length(fields);

for k=1:length(fields)
    subplot(1,nsp,k)
    nh=[healthy.(fields{k})];
    ns=[stroke.(fields{k})];
    plot([ones(length(nh),1); 2*ones(length(ns),1)],100*[nh ns]','.')
    xlim([.5 2.5])
    set(gca,'xtick',[1 2])
    set(gca,'xticklabels',{'Healthy','Stroke'})
    title(regexprep(fields(k),'_',' '))
    if k==1
        ylabel('RMS Percent Error')
    end
end

