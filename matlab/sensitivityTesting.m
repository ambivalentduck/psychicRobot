clc
clear all

%% Step 0: Set up a basic kicked movement and forward simulate.

t=0:.001:2;
coeff=calcminjerk(0,.15,0,0,0,0,0,.7);
[x,v,a]=minjerk(coeff,t);
f=zeros(length(t),2);

%forward simulate
%extract and show high accuracy.


%% Step 1: Sensitivity to gross and measured parameter misestimation.



% l1, l2, mass, x0

% whitened forces

% force sensor drift

% position noise (white)

% position drift 