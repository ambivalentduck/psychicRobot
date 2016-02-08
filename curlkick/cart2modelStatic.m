function [TBMP,TSP,TNP]=cart2modelStatic(X,Y,l1,l2,x0,mass)

% TBMP=Body Mass Proportional Torque
% TSP=Stiffness Proportional Torque
% TNP=Torque Proportional to Neither of the above

TBMP=zeros(size(X,1),2);
TSP=zeros(size(X,1),2);
TNP=zeros(size(X,1),2);

kp0=[10.8 2.83; 2.51 8.67];
kp1=[3.18 2.15; 2.34 6.18];

for k=1:size(X,1)
    %Compute real q and tau
    q=ikinStatic(X(k,1:2),l1,l2,x0);
    fJq=fJstatic(q,l1,l2);
    qdot=fJq\X(k,3:4)';
    qddot=getAlphaStatic(q,qdot,X(k,5:6)',l1,l2);
    torque=-fJq'*X(k,7:8)'; %Negative sign is consistent with extract.m
    qreal=[q; qdot; qddot];
    [D_real,C_real]=computeDCStatic(qreal(1:2),qreal(3:4),l1,l2,mass);

    %Compute desired q and tau
    q=ikinStatic(Y(k,1:2),l1,l2,x0);
    fJq=fJstatic(q,l1,l2);
    qdot=fJq\Y(k,3:4)';
    qddot=getAlphaStatic(q,qdot,Y(k,5:6)',l1,l2);
    qdes=[q; qdot; qddot];
    [D_des,C_des]=computeDCStatic(qdes(1:2),qdes(3:4),l1,l2,mass);
    
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