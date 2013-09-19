function y=extract(t,xvaf,armdynamics,x0)

global measuredVals measuredTime l1 l2

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

if nargin<4
    q0=measuredVals(1,1:4);
else
    q0=ikin(x0(1:2));
    q0(3:4)=fJ(q0)\x0(3:4)';
end

if strcmp(armdynamics,'reflex')
    [T,X]=extractionReflexHelper(t,q0);
else
    [T,X]=ode45(armdynamics,t,q0);
end
    
y=X(:,1:4);

for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2),l1,l2)*X(k,3:4)')';
end