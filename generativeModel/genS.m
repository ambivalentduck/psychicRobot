clc
clear all

M=1000;

minN=3;
meanN=5;
mu=meanN-minN;

lambda=1; %Around 10 cm/s

n=minN+floor(exprnd(mu,M,1));
figure(1)
clf
subplot(1,2,1)
[f,x]=ecdf(n);
plot(x,log(1-f))

Sn2=gamrnd(n,lambda*ones(M,1));
subplot(1,2,2)
hold on
[f,x]=ecdf(Sn2);
plot(x,f)
%hist(Sn2)
p=gamfit(Sn2)
plot(x,gamcdf(x,p(1),p(2)),'r')
a=unique(n);
for k=a'
    i=find(n==k);
    [f,x]=ecdf(Sn2(i));
    plot(x,f/(M/length(i)),'k')
end