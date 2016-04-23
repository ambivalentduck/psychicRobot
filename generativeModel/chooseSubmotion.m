function [S,xf,Uf]=chooseSubmotion(U,X,Y,U0,x0)

Q0=.001; %Cost is a loss
Q1=2*1.875^2;
alpha=.15^2*Q1; %mean peak submovement velocity = 15 cm / s

L2=(X-x0(1)).^2+(Y-x0(2)).^2;
deltaU=U-U0;
Tn2=(-deltaU-Q0)./(3*L2*Q1);
eJ=exp((-deltaU-L2.*Tn2*Q1-Q0)/alpha);
eJ(L2==0)=0; %Don't consider staying still.
eJ(Tn2<=0)=0;
eJ=eJ(:);
cumP=cumsum(eJ);
cumP=cumP/cumP(end);
chosen=find(cumP>=rand,1,'first');

if isempty(chosen)
    [~,chosen]=min(deltaU(:));
    Tn2=-Tn2;
end

S=Tn2(chosen)^-.5;
xf=[X(chosen) Y(chosen)];
Uf=U(chosen);
