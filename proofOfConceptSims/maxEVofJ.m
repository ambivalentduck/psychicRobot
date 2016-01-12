clc
clear all

global Tn2 L2 deltaU

N=1000;
L2=1;
deltaU=10;

Tn2=exprnd(1,N,1);

p=fmincon(@jobj,[0; 1],-eye(2),[0; -.3])

figure(1)
clf
hold on

Tcrit=sqrt(3*L2/deltaU);
T=linspace(sqrt(L2/deltaU),10*Tcrit);
Jdot=deltaU./T-L2.*T.^-3;
%0=deltaU-L2.*T^-2

plot(T.^-2,Jdot,'b-')
plot(deltaU/(3*L2),1,'rx')

T=(p(2)*Tn2+p(1)).^-.5;
J=deltaU./T-L2.*T.^-3;
plot(T.^-2,J,'r.','markersize',1)
