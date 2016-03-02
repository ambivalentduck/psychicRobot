function [rmseR,rmseT]=quickfitTsArray(reachStruct)

M=length(reachStruct);

N=20;
Tx=zeros(M,N);
Rx=zeros(M,N);
RoTx=zeros(M,N);

for t=1:M
    [Rx(t,:),Tx(t,:),RoTx(t,:)]=getTsArray(reachStruct(t).x,reachStruct(t).t,reachStruct(t).x0,reachStruct(t).xf,N);
end

    [~,~,~,rmseR]=fitShiftedGam(Rx(:,k));
    [~,~,~,rmseT]=fitShiftedGam(Tx(:,k).^-2);
end
