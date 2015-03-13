clc
clear all

load('free_exp_05stroke.mat')

en=dot(v',v');

figure(75)
clf
sx1=size(x,1);
ds=1:150:sx1;
plot3(x(ds,1),x(ds,2),en(ds),'.')