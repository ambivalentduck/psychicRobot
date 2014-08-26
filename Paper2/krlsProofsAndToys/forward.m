function xva=forward(t,Xd,F,K)

global measuredVals measuredTime

measuredVals=[Xd,F,K];
measuredTime=t;

[T,X]=ode45(@toyDyn,t,measuredVals(1,1:4));

ah=zeros(length(T),2);
for k=1:length(T)
    blah=toyDyn(T(k),X(k,:)')';
    ah(k,:)=blah(3:4);
end

xva=[X ah];
