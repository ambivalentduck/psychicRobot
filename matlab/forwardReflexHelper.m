function [T,Y]=forwardReflexHelper(T,Q0)

global measuredVals measuredTime errorVals errorTime inertialTorque

% Splitting samples in half allows a muscle torque lag for determining Kp that theoretically
% introduces negligable new error and/or uncertainty while still being efficient and allowing ode45.

tsim=zeros(2*length(T)-1,1);
for k=1:length(T)-1
    tsim([2*k-1 2*k])=[T(k); (T(k)+T(k+1))/2];
end
tsim(end)=T(end);

inertialTorque=0;
t=T(1);
X=Q0;

Y=zeros(length(T),6);
Y(1,:)=measuredVals(1,1:6);

for k=2:length(tsim)
    f=find((T>=(tsim(k)-.07))&(T<tsim(k))); %mildly larger buffer than necessary
    errorTime=T(f); 
    if k<4
        errorTime=0;
        errorVals=[0 0 0 0];
    else
        tNN=twoNearestNeighbor(measuredVals,measuredTime,errorTime);
        errorVals=tNN(:,1:4)-Y(f,1:4);
    end
    
    [t,X]=ode45(@armdynamicsBurdetReflexes,tsim(k-1:k),X(end,:));
    [dX,inertialTorque]=armdynamicsBurdetReflexes(t(end),X(end,:)');
    
    if mod(k,2)==1
        Y((k+1)/2,:)=[X(end,:) dX(3:4)'];
    end
end