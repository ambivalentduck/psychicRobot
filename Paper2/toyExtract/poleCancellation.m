clc
clear all

global measuredTime measuredVals M B

syms s
syms t real

tf=.6;
M=.41;
B=2.3;
K=16;
syms epsilon real
%epsilon=2;
Ke=K+epsilon;
%K=laplace(15+3*sin(10*pi*t));

ta=t/tf;
yd=.1-.1*(10*ta^3-15*ta^4+6*ta^5); %Choose a 10 cm reach
Yd=int(exp(-s*t)*yd,0,tf)

Yd=simple(Yd);

F=1/s-exp(-.1*s)/s; %A pulse from 0 to 100 ms



X=Yd-F/(M*s^2+B*s+K);
X=simple(X)

Ye=X+F/(M*s^2+B*s+Ke);

RES=(Ye-X)/(M*s^2+B*s+Ke)
checkRES=F/((M*s^2+B*s+Ke)^2)
shouldbezero=RES-checkRES
res=ilaplace(RES)
resplot=inline(vectorize(res))

figure(7)
clf
hold on
tplot=0:.001:1;
xplot=inline(vectorize(ilaplace(X)))
ydplot=inline(vectorize(yd));
F=(tplot<.1);
plot(tplot,F/20,'color',[.5 .5 .5]);
plot(tplot,xplot(tplot),'k')
quickfix=ydplot(tplot);
quickfix(tplot>tf)=0;
plot(tplot,quickfix,'g')

measuredTime=tplot;
measuredVals=zeros(length(tplot),10);
measuredVals(:,7)=1*(tplot<.1);
epsResid=.3;
measuredVals(:,9)=K+epsResid;

[T,X]=ode45(@integrateForce,tplot,zeros(4,1));
figure(89)
clf
hold on
plot(tplot,resplot(epsResid,tplot),'r')
plot(tplot,X(:,1),'b')

figure(8)
clf
hold on
epsvec=-4:.1:4;
storeme=zeros(length(epsvec(epsvec~=0)),1001);
q=0;
for eps=epsvec
    rp=resplot(eps,tplot);
    if eps~=0
        q=q+1;
        storeme(q,:)=rp;
    end
    
    plot3(tplot,eps+0*tplot,real(rp),'c')
end

figure(9)
clf
hold on
e=epsvec(epsvec~=0);
for k=1:10:length(tplot)
    plot3(e,0*e+k,real(storeme(:,k)))
end

%And the real question: is stiffness *overestimation* roughly linear?
r2=zeros(length(tplot),1);
m=r2;

reconstruct=zeros(length(tplot),1);

fi=find(e>0);
for k=1:length(tplot)
    r2(k)=corr(e(fi)',real(storeme(fi,k))).^2;
    m=[e(fi)' 0*e(fi)'+1]\real(storeme(fi,k));
    slope(k)=m(1);
    intercept(k)=m(2);
end

figure(10)
clf
subplot(2,1,1)
plot(tplot,r2,'.')
ylabel('R^2 for regression of extraction error onto stiffness overestimation')
subplot(2,1,2)
plot(tplot,slope,'.')
ylabel('Slope of that regression')
xlabel('Time since onset of movement')

figure(12)


plot(tplot,reconstruct,'.')