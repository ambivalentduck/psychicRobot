clc
clear all

load wei_2010_2.mat

for k=1:length(Subj)
    figure(k)
    clf
    [~,~,~,rmserw(k)]=fitShiftedGam(Subj(k).rw,1)
    [~,~,~,rmserts2(k)]=fitShiftedGam(Subj(k).ts2,1)
end

figure(2000)
clf
hold on
handle1=plot(ones(length(Subj),1)-.1,rmserw,'r.')
handle2=plot(ones(length(Subj),1)+.1,rmserts2,'b.')

load /home/justin/Dropbox/JustinCollabPaper/Thoroughman_2005.mat

for k=1:length(Subj)
    figure(k)
    clf
    [~,~,~,rmseTrw(k)]=fitShiftedGam(Subj(k).rw,1)
    [~,~,~,rmserTts2(k)]=fitShiftedGam(Subj(k).ts2,1)
end

figure(2000)
plot(1+ones(length(Subj),1)-.1,rmseTrw,'r.')
plot(1+ones(length(Subj),1)+.1,rmserTts2,'b.')


for k=1:8
    figure(k)
    clf
    reachStruct=processPulse(k);
    [rmseR(k),rmseT(k)]=quickfitTsArray(reachStruct)
end

figure(2000)
plot(2+ones(length(Subj),1)-.1,rmseTrw,'r.')
plot(2+ones(length(Subj),1)+.1,rmserTts2,'b.')

set(gca,'xtick',[1 2],'xticklabels',{'Wei 2010','Thoroughman 2005'})
legend([handle1, handle2],{'Rectified Work','S^{ -2}'})