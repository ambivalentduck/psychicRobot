clc
clear all

invlambdas=[1 2 4 8 16 32];
lambdas=1./invlambdas;

N=1000000;

ISIs=zeros(N,length(lambdas));

for k=1:length(lambdas)
    ISIs(:,k)=exprnd(lambdas(k),N,1);
end

[winner,column]=min(ISIs,[],2);

ensemblelambda=expfit(winner);
1/ensemblelambda
sum(invlambdas)

pwinner=zeros(1,length(lambdas));
for k=1:length(lambdas)
    pwinner(k)=sum(column==k)/length(column);
end

pwinner
pwinner-invlambdas/sum(invlambdas)