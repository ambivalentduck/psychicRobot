clc
clear all

%% Setup Cartesian space
x=-.1:.001:.1;
y=.3:.001:.5;

[X,Y]=meshgrid(x,y);

%% Setup and plot potential field
%Uf=@(X,Y) 1000*(abs(X)-.01).*(abs(X)>.01)+1000*abs(Y-.45);
Uf=@(X,Y) 20*abs(X)+8*(Y-.45).^2;
U=Uf(X,Y);

figure(1)
clf
skip=5;
surf(X(1:skip:end,1:skip:end),Y(1:skip:end,1:skip:end),U(1:skip:end,1:skip:end))

%% Perform submotion selection using predecided U and X

figure(2)
clf
hold on

M=200;
nMax=30;
ns=zeros(M,1);


for m=1:M
    C=zeros(nMax,1);
    S=C;
    L=zeros(nMax,2);
    x=zeros(nMax+1,2);
    Us=x;
    
    x(1,:)=[0,.3]+randn(1,2)*.002;
    Us(1,:)=Uf(x(1,1),x(1,2));
    n=0;
    while (norm(x(n+1,:)-[0,.45])>.005)&(n<nMax)
        n=n+1;
        [S(n),x(n+1,1:2),Us(n+1)]=chooseSubmotion(U,X,Y,Us(n),x(n,:));
        
        if n==1
            C(1)=S(n)/2;
        else
            C(n)=C(n-1)+sqrt(exprnd(.1214));
        end
        L(n,:)=x(n+1,1:2)-x(n,1:2);
    end
    C=C(1:n);
    S=S(1:n);
    L=L(1:n,:);
    t=0:.01:min(5,max(C+S));
    ydot=zeros(length(t),2);
    y=ones(length(t),1)*x(1,:);
    for k=1:n
        tau=(t'-C(k))/S(k)+.5;
        tau=max(min(tau,1),0);
        kappa=30*tau.^2-60*tau.^3+30*tau.^4;
        ydot=ydot+kappa*(L(k,:)/S(k));
        kappaPos=10*tau.^3-15*tau.^4+6*tau.^5;
        y=y+kappaPos*L(k,:);
    end
    
    plot(y(:,1),y(:,2))
    plot(y(end,1),y(end,2),'rx')
    ns(m)=n;
    %Need to calculate actual error metrics here.
end
axis equal

figure(3)
clf
hist(ns,1:25)

return

figure(4)
clf
hold on
[f,x]=ecdf(ts.^-2);
plot(x,f)
p=gamfit(ts.^-2);
plot(x,gamcdf(x,p(1),p(2)),'r')