function y=extract(t,xvaf,params,x0)

global measuredVals measuredTime fJ getAlpha

set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

measuredVals=xvaf;
measuredTime=t;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)');
    torque=fJq'*xvaf(k,7:8)';
    measuredVals(k,:)=[q' qdot' qddot' torque'];
end

if nargin<=3
    [T,X]=ode45(@armdynamics_inverted,t,measuredVals(1,1:4));
else
    q=ikin(x0(1:2));
    q0=[q fJ(q)\x0(3:4)'];
    [T,X]=ode45(@armdynamics_inverted,t,q0);
end
    
y=X(:,1:4);

for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2))*X(k,3:4)')';
end