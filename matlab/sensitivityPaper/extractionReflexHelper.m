function [T,X]=extractionReflexHelper(T,Q0)

global measuredVals measuredTime errorVals errorTime

%Deal with the first 60ms where there's presumably zero error in the
%history to cause reflexes.
f=find(T<(T(1)+.06));
if length(f)<=1
    error('You passed a single point or less?')
end
[t,X]=ode45(@armdynamicsInvertedBurdet,T(f),Q0');
sols(1).t=t';
sols(1).X=X';

lastend=f(end);
f=find((T>T(lastend))&(T<T(lastend)+.06));
while ~isempty(f)
    errorTime=sols(end).t;
    tNN=twoNearestNeighbor(measuredVals,measuredTime,errorTime);
    errorVals=(tNN(:,1:4)-sols(end).X');
    [t,X]=ode45(@armdynamicsInvertedBurdetReflexes,T([lastend f]),X(end,:));
    if length(f)~=1
        sols(end+1).t=t(2:end)';
        sols(end).X=X(2:end,:)';
    else
        sols(end+1).t=t(end)';
        sols(end).X=X(end,:)';
    end
    lastend=f(end);
    f=find((T>T(lastend))&(T<T(lastend)+.06));
end

t=[sols.t];
X=[sols.X]';