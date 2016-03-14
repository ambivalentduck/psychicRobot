clc
clear all

x=-.1:.001:.1;
y=.3:.001:.5;

[X,Y]=meshgrid(x,y);

%Uf=@(X,Y) 1000*(abs(X)-.01).*(abs(X)>.01)+1000*abs(Y-.45);
Uf=@(X,Y) sqrt(20*X.^2+8*(Y-.45).^2)
U=Uf(X,Y);

figure(1)
clf
skip=5;
surf(X(1:skip:end,1:skip:end),Y(1:skip:end,1:skip:end),U(1:skip:end,1:skip:end))

figure(2)
clf
hold on

ns=zeros(50,1);
ts=ns;

for N=1:100
    t=0;
    x=[0,.3]+randn(1,2)*.005;
    Us=Uf(x(1),x(2));
    while (norm(x(end,:)-[0,.45])>.01)&(length(t)<25)
        [t(end+1),x(end+1,1:2),Us(end+1)]=chooseSubmotion(U,X,Y,Us(end),t(end),x(end,:));
    end
    plot(x(:,1),x(:,2))
    plot(x(end,1),x(end,2),'rx')
    ns(N)=length(t);
    ts(N)=t(end);
end
axis equal

figure(3)
clf
hist(ns,1:15)

figure(4)
clf
hist(ts.^-2,20)