function [T,X]=extractionReflexHelper(T,Q0,E0)

global measuredVals measuredTime errorVals errorTime

errorVals=E0';
errorTime=T(1);

k=0;
X0=Q0;
lastf=2;
while (T(1)+.06*(k-1))<T(end)
    k=k+1;
    f=find(T<(T(1)+.06*k),1,'last');
    [t,X]=ode45(@armdynamicsInvertedBurdetReflexes,T(lastf-1:f),[X0(1:2),X0(3:4)]);
    if f-lastf>1
        sols(k).X=X(2:end-1,:)';
        errorTime=T(lastf+1:f);
    else
        sols(k).X=X(end,:)';
        errorTime=T(f);
    end
    tNN=twoNearestNeighbor(measuredVals,measuredTime,errorTime);
    errorVals=(tNN(:,1:4)-sols(k).X');
    X0=X(end,:);
    lastf=f;
end

try
    X=[Q0' sols.X]';
catch
    X=[Q0 sols.X]';
end