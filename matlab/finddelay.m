clc
clear all

figure(46)
clf
hold on

fit=zeros(10,2,8);

for S=1:8 %why not
    
    prefix='pulse';
    
    load(['../Data/Data_pulse/',prefix,num2str(S),'.mat'])
    
    a=vertcat(trials(3:479).a);
    f=vertcat(trials(3:479).f);
    
    l=size(a,1);
    
    subplot(4,2,S)
    hold on
    
    for blah=3:12
        ind=floor(blah*l/16):floor((blah+1)*l/16);
        pairs=zeros(31,2);
        for k=-30:30
            mdl=fitlm(a(ind,1),f(ind+k,1));
            pairs(31+k,:)=[k*.005 mdl.Rsquared.Ordinary];
        end
        [v,i]=max(pairs(:,2));
        plot(1000*pairs(:,1),pairs(:,2))
        plot(1000*pairs(i,1),pairs(i,2),'rx')
        fit(blah-2,:,S)=[1000*pairs(i,1),pairs(i,2)];
    end
    ave=dot(fit(:,1,S),fit(:,2,S))/sum(fit(:,2,S));
    yl=ylim;
    plot(ave*[1 1],[0 yl(2)],'m-')
    xlabel('Force Sensor Lag, ms')
    ylabel('R^2, regression of acceleration and force')
    
end


