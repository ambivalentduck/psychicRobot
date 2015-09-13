function [shift,n,T,cost]=fitTsArray(reachStruct)

M=length(reachStruct);

N=20;
Tx=zeros(M,N);
Rx=zeros(M,N);
RoTx=zeros(M,N);

for t=1:M
    [Tx(t,:),Rx(t,:),RoTx(t,:),offvec]=getTsArray(reachStruct(t).x,reachStruct(t).t,reachStruct(t).x0,reachStruct(t).xf,N);
end

oVec=mean([offvec(1:end-1) offvec(2:end)],2);
Txd=diff(Tx,1,2);
Rxd=diff(Rx,1,2);

shift=zeros(N-1,1);
n=zeros(N-1,1);
T=zeros(N-1,1);
shiftCum=zeros(N,1);
nCum=zeros(N,1);
TCum=zeros(N,1);
for k=1:N-1
    [shift(k),n(k),T(k)]=fitShiftedGam(Rxd(:,k).*(Txd(:,k).^-2));
    [shiftCum(k),nCum(k),TCum(k)]=fitShiftedGam(Tx(:,k).^-2);
end
k=N;
figure(2)
clf
Tx(:,k)
[shiftCum(k),nCum(k),TCum(k)]=fitShiftedGam(Tx(:,k).^-2,1);


figure(1)
clf
subplot(1,3,1)
hold on
plot(oVec,shift,'b-')
plot(offvec,shiftCum,'g-')
ylabel('Param value')
subplot(1,3,2)
hold on
plot(oVec,n)
plot(offvec,nCum,'g-')
subplot(1,3,3)
hold on
plot(oVec,T)
plot(offvec,TCum,'g-')

return
subplot(2,3,4)
loglog(offvec,gradient(shift)./gradient(offvec))
xlabel('Shift')
ylabel('Spatial Gradient')
subplot(2,3,5)
plot(offvec,gradient(n)./gradient(offvec))
xlabel('n')
subplot(2,3,6)
plot(offvec,gradient(T)./gradient(offvec))
xlabel('T')


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

