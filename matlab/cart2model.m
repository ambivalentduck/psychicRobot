function [SNE,OUT,BOUT]=cart2model(X,Y)

% Assume that forces perpendicular to the direction of movement don't
% *quickly* effect parallel progress
% (X is x,y and derivatives measured; Y is y-component desired/typical)

% SNE=sum of non-error-dependent terms
% Ep=position error
% Ev=velocity error
% Tm=muscle torque magnitude vector, needed for Perrault-Franklin style fitting
% BOUT=[Kp0*(Ep+Ev/12) Kp1*(Tm.*(Ep+Ev/12))]
% OUT=[Ep Ev Ep.*Tm Ev.*Tm]
% In theory, Kp_burdet = BOUT\SNE, Kp_burdet = [Kp0 Kp1]

global fJ getAlpha

SNE=zeros(2,size(X,1));
Ep=SNE;
Ev=SNE;
Tm=SNE;

kp0=[10.8 2.83; 2.51 8.67];
kp1=[3.18 2.15; 2.34 6.18];

for k=1:size(X,1)
    q=ikin(X(k,1:2));
    fJq=fJ(q);
    qdot=fJq\X(k,3:4)';
    qddot=getAlpha(q,qdot,X(k,5:6)');
    torque=-fJq'*X(k,7:8)'; %This line could scrooge it... with the -sign
    qreal=[q; qdot; qddot];
    [D_real,C_real]=computeDC(qreal(1:2),qreal(3:4));
    
    q=ikin([X(k,1) Y(k,1)]);
    fJq=fJ(q);
    qdot=fJq\[X(k,3) Y(k,2)]';
    qddot=getAlpha(q,qdot,[X(k,5) Y(k,3)]');
    qdes=[q; qdot; qddot];
    [D_des,C_des]=computeDC(qdes(1:2),qdes(3:4));
    
    Tff=D_des*qdes(5:6)+C_des;
    Tinertial=D_real*qreal(5:6)+C_real;
    SNE(:,k)=Tff-Tinertial-torque;
    
    Ep(:,k)=qreal(1:2)-qdes(1:2);
    Ev(:,k)=qreal(3:4)-qdes(3:4);
    Tm(:,k)=abs(Tinertial+torque);
end

SNE=SNE';
Ep=Ep';
Ev=Ev';
Tm=Tm';

OUT=[Ep Ev Ep.*Tm Ev.*Tm];
BOUT=[(Ep+Ev/12)*kp0' (Tm.*(Ep+Ev/12))*kp1'];
[BOUT(:,[1 3]);BOUT(:,[2 4])]\[SNE(:,1); SNE(:,2)]