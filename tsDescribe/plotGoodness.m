clc
clear all

load wei_2010_2.mat

for k=1:length(Subj)
    figure(k)
    clf
    fitShiftedGam(Subj(k).rw,1)
end