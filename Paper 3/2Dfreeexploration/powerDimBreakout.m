clc
clear all
close all

k=4;
sublets={'JL','JT','LP','MF','YM'};

load(['free_exp_',sublets{k},'.mat'])

Px=a(:,1).*v(:,1);
Py=a(:,2).*v(:,2);
P=dot(a',v')';

figure(1)
hold on
ecdf(P)
ecdf(Px)
ecdf(Py)

[counts,bins]=hist([P,Px,Py],30);
figure(2)
plot(abs(bins),log(counts),'.')
