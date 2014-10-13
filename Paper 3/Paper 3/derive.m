clc
clear all

syms m x t dummy real

x=10*t^3-15*t^4+6*t^5
xdot=diff(x,t)
xddot=diff(xdot,t)

%Know U(t)
U=int(xddot*xdot,t)

ivU=inline(vectorize(U));
z=0:.01:1;
figure(1)
clf
subplot(2,2,1)
plot(z,-ivU(z))
title('Potential vs Time')
ylabel('Virtual Potential, Energy Units')
xlabel('Time')

%Know x(t)
xt = @(t) 6*t.^5-15*t.^4+10*t.^3;
subplot(2,2,2)
plot(z,xt(z))
title('1D Position vs Time')
ylabel('Position, Distance Units')
xlabel('Time')

subplot(2,2,3:4)
X=xt(z);
P=-ivU(z);
hold on
plot(X,P)

for k=2
    p=polyfit(X,P,k);
    plot(X,polyval(p,X),'k.')
end
legend('Exact','4th Order Polynomial Fit')
title('Virtual Potential vs State')
ylabel('Virtual Potential, Energy Units')
xlabel('Position, Distance Units')