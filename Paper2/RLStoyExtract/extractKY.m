function [xva,tsolved,Ksolved]=reverse(t,Xh,F,K0)

global measuredVals measuredTime

measuredVals=[Xh,F,K0+0*t']; %Have to seed K
measuredTime=t;

Ksolved=t;
Ksolved(1)=K0;
nT=5;
nK=nT*2; %We want an overdefined matrix
maxKRatio=.3;

for k=1:length(t)-nT %Demand that 
    tspan=t(k):t(k+nT)
    minK=Ksolved(k);
    for m=1:nK
        measuredVals(k:end,9)=(1+maxKRatio*(m-1)/(nK-1))*minK;
        [T,X]=ode45(@toyInvDyn,tspan,measuredVals(1,1:4));

ad=zeros(length(T),2);
for k=1:length(T)
    blah=toyInvDyn(T(k),X(k,:)')';
    ad(k,:)=blah(3:4);
end

xva=[X ad];

tsolved=t(1:end-5