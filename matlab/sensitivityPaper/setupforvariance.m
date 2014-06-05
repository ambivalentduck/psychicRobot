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

for k=[] %6:10
    blah=find(t>.7);
    w=randn(2*length(t),2);
    [Fb,Fa]=butter(4,5*2*pi*.005,'low'); %6-10
    %[Fb,Fa]=butter(10,3*.005,'low'); %1-5
    f=2.7*1.5*filter(Fb,Fa,w);
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

f=zeros(length(t),2);
k=10;
for force=[-15 15]
    for onset=([.1 .5]*.15)
        k=k+1;
        first=find(x(:,1)>=onset,1,'first');
        last=find(t>=(t(first)+.15),1,'first');
        fspecific=f;
        fspecific(first:last,2)=force;

        xvaf=[x v a fspecific];
        xsim=forwardSim(paramsPopulator,t,xvaf);

        figure(k)
        clf
        hold on
        plot(x(:,1),x(:,2),'b',xsim(:,1),xsim(:,2),'k.')
        quiver(xsim(:,1),xsim(:,2),fspecific(:,1),fspecific(:,2),'b')
        axis equal

        yex=extract(t,[xsim fspecific],'reflex');
        plot(yex(:,1),yex(:,2),'r.')
        plot(x(:,1),x(:,2),'b')

        legend('Intent','Forward Sim','Forces','Extracted Intent')

        vals(k).t=t;
        vals(k).xvaf=[xsim fspecific];
    end
end

for k=11:14
    simsforvariance(num2str(k),vals(k).t,vals(k).xvaf)
end