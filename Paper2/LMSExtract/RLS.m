function [w,p]=RLS(x,y,wm1,pm1,lambda)

k=(lambda*pm1*x)/(1+lambda*x'*pm1*x);
alpha=y-x'*wm1;

w=wm1+k*alpha;
p=lambda*pm1-lambda*k*x'*pm1;