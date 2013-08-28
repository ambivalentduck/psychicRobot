function y=extract(t,xvaf,params,armdynamics)

global measuredVals measuredTime fJ getAlpha

set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

measuredVals=xvaf;
measuredTime=t;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)');
    torque=-fJq'*xvaf(k,7:8)';
    measuredVals(k,:)=[q' qdot' qddot' torque'];
end

if strcmp(armdynamics,'reflex')
    [T,X]=extractionReflexHelper(t,measuredVals(1,1:4));
else
    [T,X]=ode45(armdynamics,t,measuredVals(1,1:4));
end
    
y=X(:,1:4);
size(y)
length(T)

for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2))*X(k,3:4)')';
end