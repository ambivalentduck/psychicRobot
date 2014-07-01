clc
clear all

load OATdat_all.mat

t=0:.005:2;
coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
tcalc=t;
tcalc(t>=.7)=.7;
[x,v,a,j]=minjerk(coeff,tcalc);
x=x';
v=v';
a=a';
j=j';

% f=zeros(length(t),2);
% 
% f((t>=.1)&(t<=.15),2)=15;
f=white.xvaf(:,7:8);

xvaf=[x v a f];

bins=0:.005:.15;
[trash,ref]=getMUE(bins,0*bins,xvaf);

in_for=xvaf;
xsim=forwardSim(paramsPopulator,t,in_for);
jsim=[gradient(xsim(:,5))/(t(2)-t(1)) gradient(xsim(:,6))/(t(2)-t(1))];

figure(1)
clf
hold on
plot(x(:,1),x(:,2),'b',xsim(:,1),xsim(:,2),'k.')
quiver(xsim(:,1),xsim(:,2),f(:,1),f(:,2),'b')
axis equal

freqs=logspace(-2,2,30);
for k=1:length(freqs)
    k
    yex(k).val1=extract(t,[xsim f(:,1) f(:,2)+cos(2*pi*freqs(k)*t')],'reflex');
    yex(k).val2=extract(t,[xsim f(:,1) f(:,2)+2*cos(2*pi*freqs(k)*t')],'reflex');
    plot(yex(k).val1(:,1),yex(k).val1(:,2),'r.')
    yex(k).mue1=getMUE(bins,ref,yex(k).val1);
    yex(k).mue2=getMUE(bins,ref,yex(k).val2);
    drawnow
end

figure(2)
clf
loglog(freqs,100*[yex.mue1],'b',freqs,100*[yex.mue2],'r')
xlabel('Error frequency, Hz')
ylabel('MUE, cm')

figure(3)
clf
hold on
plot(x(:,1),x(:,2),'b',xsim(:,1),xsim(:,2),'k.')
quiver(xsim(:,1),xsim(:,2),f(:,1),f(:,2),'b')
axis equal

freqs=logspace(-4,2,20);
for k=1:length(freqs)
    k
    yex(k).valj=extract(t,[xsim f+freqs(k)*jerror],'reflex');
    plot(yex(k).valj(:,1),yex(k).valj(:,2),'r.')
    yex(k).muej=getMUE(bins,ref,yex(k).valj);
    drawnow
end

yex=yex(1:20)

figure(4)
clf
loglog(freqs,100*[yex.muej])
ylabel('MUE, cm')
xlabel('Jerk-order Impedance, N s^3 m^{-1}','interpreter','tex')
