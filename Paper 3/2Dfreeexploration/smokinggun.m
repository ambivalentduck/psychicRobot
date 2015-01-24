clc
clear all

load('free_exp_05stroke.mat')

en=dot(v',v');

figure(75)
clf
plot3(x(:,1),x(:,2),en,'.')