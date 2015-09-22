clc
clear all

figure(1)
clf
syms t y real

tc=.5;
ts=1;
ta=(t-tc)/ts+.5;
tK=0:.001:1;
color='mgrc';

subplot(3,1,1)
hold on
K=(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
s=solve(K-y,t)
fk=inline(vectorize(K));
plot(tK,fk(tK),'b')

yv=0:.001:1.875;
for k=1:4;
    fs=inline(vectorize(s(k)))
    plot(fs(yv),yv,'color',color(k))
end
title('Velocity')
pretty(s(1))


subplot(3,1,2)
hold on
K=.5*((30*ta.^2-60*ta.^3+30*ta.^4)/ts)^2;
s=solve(K-y,t)
fk=inline(vectorize(K))
plot(tK,fk(tK),'b')
yE=0:.001:(.5*1.875^2);
for k=1:4;
    fs=inline(vectorize(s(k)))
    plot(fs(yE),yE,'color',color(k))
end
pretty(s(1))

syms t y positive
subplot(3,1,3)
hold on
K=((60*ta-180*ta.^2+120*ta.^3).*(30*ta.^2-60*ta.^3+30*ta.^4))/ts^3;
K=simplify(K);
dK=diff(K,t);
solve(dK)
fk=inline(vectorize(K))

assumeAlso(y<fk(1/2 - 7^(1/2)/14))
assumeAlso(y>0)
assumeAlso(t>0)
assumeAlso(t<1/2 - 7^(1/2)/14)

s=solve(K-y,t,'MaxDegree',4)

plot(tK,fk(tK),'b')
yP=0:.001:(.5*1.875^2);
for k=1:4;
    fs=inline(vectorize(s(k)))
    plot(fs(yP),yP,'color',color(k))
end
pretty(s(1))