function [shift,n,A]=fitTsArray(reachStruct)

M=length(reachStruct);

N=20;
Tx=zeros(M,N);
Rx=zeros(M,N);
RoTx=zeros(M,N);

for t=1:M
    [Rx(t,:),Tx(t,:),RoTx(t,:),offvec]=getTsArray(reachStruct(t).x,reachStruct(t).t,reachStruct(t).x0,reachStruct(t).xf,N);
end

oVec=mean([offvec(1:end-1) offvec(2:end)],2);
Txd=diff(Tx,1,2);
Rxd=diff(Rx,1,2);

shift=zeros(N-1,1);
n=zeros(N-1,1);
A=zeros(N-1,1);
shiftCumR=zeros(N,1);
nCumR=zeros(N,1);
ACumR=zeros(N,1);

shiftCumT=zeros(N,1);
nCumT=zeros(N,1);
ACumT=zeros(N,1);

for k=1:N-1
    %[shift(k),n(k),A(k)]=fitShiftedGam(Rxd(:,k));
    [shiftCumR(k),nCumR(k),ACumR(k)]=fitShiftedGam(Rx(:,k));
    [shiftCumT(k),nCumT(k),ACumT(k)]=fitShiftedGam(Tx(:,k).^-2);
end
k=N;
figure(2)
clf
subplot(2,1,1)
[shiftCumR(k),nCumR(k),TCumR(k)]=fitShiftedGam(Rx(:,k),1);
subplot(2,1,2)
[shiftCumT(k),nCumT(k),ACumT(k)]=fitShiftedGam(Tx(:,k).^-2,1);

k=5;
figure(5)
clf
subplot(2,1,1)
[shiftCumR(k),nCumR(k),TCumR(k)]=fitShiftedGam(Rx(:,k),1);
subplot(2,1,2)
[shiftCumT(k),nCumT(k),ACumT(k)]=fitShiftedGam(Tx(:,k).^-2,1);


figure(1)
clf
Ny=2;
Nx=3;
subplot(Ny,Nx,1)
hold on
plot(offvec,shiftCumR,'b-')
ylabel('Param value')
subplot(Ny,Nx,2)
hold on
plot(offvec,nCumR,'b-')
subplot(Ny,Nx,3)
hold on
plot(offvec,ACumR,'b-')
subplot(Ny,Nx,4)
hold on
plot(offvec,shiftCumT,'g-')
ylabel('Param value')
subplot(Ny,Nx,5)
hold on
plot(offvec,nCumT,'g-')
subplot(Ny,Nx,6)
hold on
plot(offvec,ACumT,'g-')

return


figure(SUB)
clf
subplot(1,2,1)
hold on
tessesfixed=abs(tesses(tesses~=-1));
[shift,n,T]=fitShiftedGam(tessesfixed);
x=tessesfixed.^-2;
cost=gamlike([n T],x-shift);
x=x(x<35); %Not censoring, just visibility for plotting

nbins=25;
edges=linspace(min(x),max(x),nbins+1)';
counts=zeros(nbins,1);
for k=2:length(edges)
    counts(k-1)=sum((x>=edges(k-1))&(x<=edges(k)));
end
counts=counts/sum(counts);

centers=(edges(2:end)+edges(1:end-1))/2;
bar(centers,counts)
gcx=gamcdf(edges-shift,n,T);
plot(centers,gcx(2:end)-gcx(1:end-1),'r')
ylabel('Normalized Frequency')
xlabel('t_s^{-2}')
title('Histogram of t_s^{-2}')
subplot(1,2,2)
hold on
ecdf(x,'bounds','on')
sortX=sort(x);
plot(sortX,gamcdf(sortX-shift,n,T),'r')
title(['U=',num2str(shift),' n=',num2str(n),' T=',num2str(T)])

