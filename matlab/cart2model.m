function [TBMP,TSP,TNP]=cart2model(X,Y)

% Assume that forces perpendicular to the direction of movement don't
% *quickly* effect parallel progress
% (X is x,y and derivatives measured; Y is y-component desired/typical)

% TBMP=Body Mass Proportional Torque
% TSP=Stiffness Proportional Torque
% TNP=Torque Proportional to Neither of the above

global fJ getAlpha

TBMP=zeros(size(X,1),2);
TSP=zeros(size(X,1),2);
TNP=zeros(size(X,1),2);

kp0=[10.8 2.83; 2.51 8.67];
kp1=[3.18 2.15; 2.34 6.18];

for k=1:size(X,1)
    %Compute real q and tau
    q=ikin(X(k,1:2));
    fJq=fJ(q);
    qdot=fJq\X(k,3:4)';
    qddot=getAlpha(q,qdot,X(k,5:6)');
    torque=-fJq'*X(k,7:8)'; %Negative sign is consistent with extract.m
    qreal=[q; qdot; qddot];
    [D_real,C_real]=computeDC(qreal(1:2),qreal(3:4));

    %Compute desired q and tau
    q=ikin([X(k,1) Y(k,1)]);
    fJq=fJ(q);
    qdot=fJq\[X(k,3) Y(k,2)]';
    qddot=getAlpha(q,qdot,[X(k,5) Y(k,3)]');
    qdes=[q; qdot; qddot];
    [D_des,C_des]=computeDC(qdes(1:2),qdes(3:4));
    
    %Find bodymass-proportional, stiffness-proportional, and neither-proportional terms
    Tff=D_des*qdes(5:6)+C_des;
    Tinertial=D_real*qreal(5:6)+C_real;
    Tmuscle=abs(Tinertial+torque);
    Ep=qreal(1:2)-qdes(1:2);
    Ev=qreal(3:4)-qdes(3:4);
    
    %TNP=m(TBMP)+kTSP
    TNP(k,:)=torque';
    TBMP(k,:)=(Tff-Tinertial)';
    TSP(k,:)=-((kp0+kp1*diag(Tmuscle))*(Ep+Ev/12))';
end