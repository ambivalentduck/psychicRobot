clc
clear all

t=0:.005:2;
coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
tcalc=t;
tcalc(t>=.7)=.7;
[x,v,a]=minjerk(coeff,tcalc);
x=x';
v=v';
a=a';

for k=1:5
    blah=find(t>.7);
    w=randn(2*length(t),2);
    [Fb,Fa]=butter(10,3*.005,'low');
    f=7*filter(Fb,Fa,w);
    f=f(length(t)+1:end,:);
    f(blah:end,:)=0; %turn it off when you get close

    xvaf=[x v a f];

    measuredVals=xvaf;
    measuredTime=t;

    xsim=forwardSim(paramsPopulator,t,xvaf);

    figure(k)
    clf
    hold on
    plot(x(:,1),x(:,2),'b',xsim(:,1),xsim(:,2),'k.')
    quiver(xsim(:,1),xsim(:,2),f(:,1),f(:,2),'b')
    axis equal

    yex=extract(t,[xsim f],'reflex');
    plot(yex(:,1),yex(:,2),'r.')
    plot(x(:,1),x(:,2),'b')

    legend('Intent','Forward Sim','Forces','Extracted Intent')

    xvaf=[xsim f];
    vals(k).t=t;
    vals(k).xvaf=xvaf;
end

for k=1:5
    simsforvariance(num2str(k),vals(k).t,vals(k).xvaf)
end