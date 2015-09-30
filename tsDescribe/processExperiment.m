clc
clear all

%load kording_2004
load wei_2010_1

figure(1)
clf
subplot(1,4,1)
boxplot(n)
ylabel('N')

subplot(1,4,2)
boxplot(shift)
ylabel('U')

subplot(1,4,3)
boxplot(T)
ylabel('A')

for k=1:length(n)
    X=ts{k}.^-2;
    r=randperm(length(X));
    train=X(r(1:floor(r/2)));
    test=X(r(floor(r/2)+1:end));
%     train=X;
%     test=X;
    [U,N,A]=fitShiftedGam(train); 
    sgc=@(x) shiftedGamCDF(x,U,N,A);
    [F,x]=ecdf(train);
    Strain(k)=std(F-sgc(x))/sqrt(length(x));
    [F,x]=ecdf(test);
    Stest(k)=std(F-sgc(x))/sqrt(length(x));
end


serror=[Strain' Stest']*100

subplot(1,4,4)
boxplot(serror,'labels',{'Train','Test'})


