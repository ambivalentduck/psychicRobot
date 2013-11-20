function y=extract2(t,xvaf,armdynamics,brkp,gains)

global measuredVals measuredTime fJ getAlpha kpgain

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

q0=measuredVals(1,1:4);

for k=2:length(brkp)
    kpgain=gains(k-1);
    tsub=t((t>=brkp(k-1))&(t<brkp(k)));
    if strcmp(armdynamics,'reflex')
        [T,x]=extractionReflexHelper(tsub,q0);
    else
        [T,x]=ode45(armdynamics,tsub,q0);
    end
    vals(k-1).x=x';
    q0=x(end,1:4)';
end

X=[vals.x]';

y=X(:,1:4);

for k=1:size(X,1)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2))*X(k,3:4)')';
end