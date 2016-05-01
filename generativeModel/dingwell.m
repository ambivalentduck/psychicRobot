clc
clear all

skews=[];

for k=1:9
    
    c=0;
    load(['../Data/Dingwell/YH0',num2str(k),'.mat'])

    for spd=1:5
        for tr=1:2
            eval(['X=SPD',num2str(spd),'TR',num2str(tr),'(:,2);']);
            c=c+1;
            skews(k,c)=skewness(X);
        end
    end
end

figure(1)
boxplot(skews')