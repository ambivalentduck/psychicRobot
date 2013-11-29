clc
clear all

% load all_white_workspace.mat
% 
% white.oat=OAT;
% white.xvaf=xvaf;
% white.yex=yex;
% 
% load baselines.mat
% load OAT_KICK.mat
% load KICK.mat
% load KICKsm.mat
% 
% kick.oat=OAT;
% 
% t=0:.005:2;
% coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
% tcalc=t;
% tcalc(t>=.7)=.7;
% [x,v,a]=minjerk(coeff,tcalc);
% x=x';
% v=v';
% a=a';
% 
% f=zeros(length(t),2);
% 
% f((t>=.1)&(t<=.15),2)=15;
% 
% xvaf=[x v a f];
% 
% kick.xvaf=xvaf;
% kick.yex=yex;
% 
% save('OATdat_all.mat','kick','white')

load('OATdat_all.mat')

