function [R2,A]=plotExp(bins,counts)

f=find(counts~=0);
bins=bins(f);
y=log(counts(f))';
X=[ones(length(bins),1) bins'];

b=X\y;
yhat=X*b;

A=b(2);
R2 = 1 - sum((y - yhat).^2)/sum((y - mean(y)).^2);

plot(bins',yhat,'r-')