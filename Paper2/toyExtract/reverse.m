function xva=reverse(t,Xh,F,K)

global measuredVals measuredTime

measuredVals=[Xh,F,K];
measuredTime=t;

[T,X]=ode45(@toyInvDyn,t,measuredVals(1,1:4));

ad=zeros(length(T),2);
for k=1:length(T)
    blah=toyInvDyn(T(k),X(k,:)')';
    ad(k,:)=blah(3:4);
end

xva=[X ad];
