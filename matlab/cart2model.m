function [SNE,C0,C1]=cart2model(X,Y)

% Assume that forces perpendicular to the direction of movement don't
% *quickly* effect parallel progress
% (X is x,y and derivatives measured; Y is y-component desired/typical)

% SNE=sum of non-error-dependent terms
% Ep=position error
% Ev=velocity error
% Tm=torque vector, needed for Perrault-Franklin style fitting
% C0=for fitting a constant gain to the Torque independent stiffness of Burdet
% C1=for fitting a constant gain to the Torque independent stiffness of Burdet


global fJ getAlpha

S=zeros(2,size(X,1));
C0=S;
C1=S;

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
    S(:,k)=Tff-Tinertial-torque;
    
    error=(qreal(1:2)-qdes(1:2))+(qreal(3:4)-qdes(3:4))/12;
    C0(:,k)=kp0*error;
    C1(:,k)=kp1*diag(abs(Tinertial+torque))*error;
end

S=S';
C0=C0';
C1=C1';