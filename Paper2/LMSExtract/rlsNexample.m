clc
clear all

N=50;

S2N=.01;

M=3;
wReal=[-1.5 2 1]
x=randn(M,N);
xt=x;
%xt(2,:)=sin(xt(2,:));
Y=wReal*xt;
intent=cos(12.7*(1:N))+5;
y=Y+S2N*std(Y)*randn(1,N)+intent;
linRegress=(x'\y')' %#ok<*NOPTS>

w=zeros(M,N);
lambda=.99999;
P=zeros(M,M,N);
P(:,:,1)=eye(M);

c=2;

while c<=N
%     k=(lambda*P(:,:,c-1)*x(:,c))/(1+lambda*x(:,c)'*P(:,:,c-1)*x(:,c));
%     alpha=y(c)-x(:,c)'*w(:,c-1);
%    
%     w(:,c)=w(:,c-1)+k*alpha;
%     P(:,:,c)=lambda*P(:,:,c-1)-lambda*k*x(:,c)'*P(:,:,c-1);
    [w(:,c),P(:,:,c)]=RLS(x(:,c),y(:,c),w(:,c-1),P(:,:,c-1),lambda);
    
    c=c+1;
end
w_end=w(:,end)'

close all
figure(1)
subplot(3,1,1)
yp=sum(w.*x);
plot(1:N,y,'r',1:N,Y,'b',1:N,yp,'g')
legend('With noise','Pure','Predicted')

subplot(3,1,2)
plot(1:N,abs(Y-yp))

subplot(3,1,3)
plot(1:N,intent,1:N,y-yp)

