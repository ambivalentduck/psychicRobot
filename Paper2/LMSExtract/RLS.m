function [w,p]=RLS(x,y,wm1,pm1,lambda)

k=(lambda*pm1*x)/(1+lambda*x'*pm1*x);
alpha=y-x'*wm1;

w=wm1+k*alpha;
p=lambda*pm1-lambda*k*x'*pm1;

%     k=(lambda*P(:,:,c-1)*x(:,c))/(1+lambda*x(:,c)'*P(:,:,c-1)*x(:,c));
%     alpha=y(c)-x(:,c)'*w(:,c-1);
%    
%     w(:,c)=w(:,c-1)+k*alpha;
%     P(:,:,c)=lambda*P(:,:,c-1)-lambda*k*x(:,c)'*P(:,:,c-1);
