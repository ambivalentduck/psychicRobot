clc
clear all

N=1000;
X=exprnd(2,N);

Y=X(:,1)-X(:,2);

figure(67)
clf
hold on
[Fx,x]=ecdf(X(:,1));
[Fy,y]=ecdf(abs(Y));
plot(x,Fx,'b',y,Fy,'r')