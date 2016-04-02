function [tf,xf,Uf]=chooseSubmotion(U,X,Y,U0,t0,x0)

B=.15^2; %mean peak submovement velocity = 15 cm / s
Q0=-5; %Cost is a loss
Q1=-2*1.875^2;

L2=(X-x0(1)).^2+(Y-x0(2)).^2;
deltaU=U-U0;
Tn2=(-deltaU+Q0)./(3*L2*Q1);
eJ=exp(-(deltaU+L2.*Tn2*Q1+Q0)/B);
eJ(L2==0)=0; %Don't consider staying still.
eJ(Tn2<=0)=0;
eJ=eJ(:);
cumP=cumsum(eJ);
cumP=cumP/cumP(end);
chosen=find(cumP>=rand,1,'first');

tf=t0+Tn2(chosen)^-.5;
xf=[X(chosen) Y(chosen)];
Uf=U(chosen);
