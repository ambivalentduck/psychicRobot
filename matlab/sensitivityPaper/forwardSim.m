function [y,ysm]=forwardSim(params,t,xvaf)

global l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kp kd kpgain kpkdratio reflexratio kpkdreflexratio measuredVals measuredTime

v=params;
l1=v(1);
l2=v(2);
lc1=l1*v(3);
lc2=l2*v(4);
mass=v(5);
m1=v(6)*mass;
m2=v(7)*mass;
I1=m1*(v(8)*l1)^2;
I2=m2*(v(9)*l2)^2;

x0=[v(10) v(11)];
xvafnom=xvaf;
xvaf(:,7:8)=[xvafnom(:,7)+v(12)+randn(length(t),1)*v(14) xvafnom(:,8)+v(13)+randn(length(t),1)*v(15)];

kpgain=v(16);
kp0=v(17)*[v(21) v(22); v(23) v(24)];
kp1=v(18)*[v(25) v(26); v(27) v(28)];
kp=v(19)*[v(29) v(30); v(31) v(32)];
kd=v(20)*[v(33) v(34); v(35) v(36)];

kpkdratio=v(37);
reflexratio=v(38);
kpkdreflexratio=v(39);

measuredVals=xvaf;
measuredTime=t;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q,l1,l2);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)',l1,l2);
    torque=-fJq'*xvaf(k,7:8)';
    measuredVals(k,:)=[q' qdot' qddot' torque'];
end

q0=measuredVals(1,1:4);

%forward simulate
[T,X]=forwardReflexHelper(t,q0);
[T,XSM]=ode45(@armdynamicsShadMuss,t,q0);

y=zeros(length(T),6);
ysm=y;

for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2),l1,l2)*X(k,3:4)')';
    ysm(k,1:2)=fkin(XSM(k,1:2));
    ysm(k,3:4)=(fJ(XSM(k,1:2),l1,l2)*XSM(k,3:4)')';
end

gT=gradient(t)';
y(:,5)=gradient(y(:,3))./gT;
y(:,6)=gradient(y(:,4))./gT;
ysm(:,5)=gradient(ysm(:,3))./gT;
ysm(:,6)=gradient(ysm(:,4))./gT;