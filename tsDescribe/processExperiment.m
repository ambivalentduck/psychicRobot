clc
clear all

%load kording_2004
load wei_2010_1

figure(1)
clf
subplot(1,3,1)
boxplot(n)

subplot(1,3,2)
boxplot(shift)

subplot(1,3,3)
boxplot(T)

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
figure(1)
clf
hold on
[F,x,Fl,Fu]=ecdf(train);
h=fill([x(~isnan(Fl)); wrev(x(~isnan(Fu)))],[Fl(~isnan(Fl)); wrev(Fu(~isnan(Fu)))],'k');
set(h,'edgealpha',0,'facealpha',.5)
plot(x,F,'k')
plot(x,sgc(x),'r')

[Strain' Stest']*100


