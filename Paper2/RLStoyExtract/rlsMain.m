clc
clear all

global M B

%% Constant K, prove tautology and other housekeeping
tspan=.5;
ta=(0:.001:tspan)/tspan;

xi=[1; 0];
xf=[0; 0];

xd=(xi*ones(size(ta))+(xf-xi)*(10*ta.^3-15*ta.^4+6*ta.^5))';
vd=((xf-xi)*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';
ad=((xf-xi)*(60*ta-180*ta.^2+120*ta.^3)/(tspan^2))';
% tshort=0:.001:tspan;
% xd=([1; 0]*sin(tshort))';
% vd=([1; 0]*cos(tshort))';
% ad=([1; 0]*-sin(tshort))';

t=0:.001:1.5;
Xd=[[xd vd ad];zeros(length(t)-length(ta),6)];
Xd(:,1)=Xd(:,1)+1;
F=zeros(length(t),2);
F((t>.1)&(t<.2),:)=-10;
K=15*ones(length(t),2);

figure(1)
clf
hold on
plot(t,Xd(:,1),'g')

M=.41;
B=2.3;

Xh=forward(t,Xd,F,K);
plot(t,Xh(:,1),'k')

Xx=reverse(t,Xh,F,K);
plot(t,Xx(:,1),'r.')

arrow([t' Xh(:,1)],[t' Xh(:,1)]+.05*[0*t' F(:,1)],.6*[1 1 1],.2);

title('Proof of Tautology Being Intact')
xlabel('Time, Seconds')
ylabel('Displacement in X')
legend('Ground Truth Intention','Hand (Stiffness=15)','Intention (Stiffness=15)','Forces')

%% Sobol spanned, short-time extraction bursts

sParams=1;
nCandidates=5;
Kml=16;
Kspan=.2;

p=sobolset(sParams,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

candidateK=Kml+Kml*Kspan*(p(1:nCandidates,:)-.5); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.

%You know you have some settling time, so figure out what it is.

figure(2)
clf
subplot(2,1,1)
hold on
plot(t,Xd(:,1),'g')

Xml=reverse(t,Xh,F,Kml+K*0);
plot(t,Xml(:,1),'r')

Xs=zeros(length(t),nCandidates);
Ys=zeros(length(t),nCandidates);
for k=1:nCandidates
    cand{k}=reverse(t,Xh,F,candidateK(k)+K*0);
    Xs(:,k)=cand{k}(:,1);
    Ys(:,k)=cand{k}(:,2);
    plot(t,Xs(:,k),'b')
end

regressX=candidateK\Xs';
regressY=candidateK\Ys';

P=eye(3);
w=zeros(3,2,length(t));

subplot(2,1,2)
hold on
for k=1:length(t)
    [w(:,:,k+1),P]=RLS([regressX(k);regressY(k);1],Xml(k,1:2)',w(:,:,k),P,.99);
    if ~mod(k,10)
    plot(t(k),w(1,1,k),'r.')
    plot(t(k),w(2,2,k),'b.')
    drawnow
    end
end

figure(3)
clf
hold on
plot(Xd(:,1),Xd(:,2),'g.')
plot(squeeze(w(3,1,:)),squeeze(w(3,2,:)),'b')


%% Break there

return

figure(2)
clf
subplot(2,1,1)
hold on
plot(t,Xd(:,1),'g')
%arrow([t' Xd(:,1)],[t' Xd(:,1)]+.05*[0*t' F(:,1)],.6*[1 1 1],.2);

Xh=forward(t,Xd,F,K);
plot(t,Xh(:,1),'k')

Klow=13;
Khigh=16;

XxL=reverse(t,Xh,F,Klow+0*K);
plot(t,XxL(:,1),'r')

XxH=reverse(t,Xh,F,Khigh+0*K);
plot(t,XxH(:,1),'r')
title('You can recover the twue intent by overestimation correction')
ylabel('Displacement in X')
xlim([0 1.5])
legend('Ground Truth Intention','Hand (Stiffness=15)','Extraction with overestimated stiffness (Stiffness=18)','Extraction with overestimated stiffness (Stiffness=19)')

subplot(2,1,2)
hold on
E=XxL(:,1)-XxH(:,1);
plot(t,E,'c')

plot(t,Xd(:,1),'g')
eps=E\XxL(:,1)
plot(t,XxL(:,1)-eps*E,'r')
legend('Difference of Overestimates','Ground Truth Intention','Corrected Extraction')
title(['Solved Stiffness Overestimation: Eps = ',num2str(-eps*(Khigh-Klow))])
xlabel('Time, Seconds')
ylabel('Displacement in X')


%% Time varying K

Ksmooth=40+5*sin(2*pi*t/.25+rand)'; %smooth(5+20*rand(length(t),1),30);
Xh=forward(t,Xd,F,[Ksmooth K(:,2)]);
Klow=9;
Kmid=12;
Khigh=12.1;

XxL=reverse(t,Xh,F,Klow+0*K);
XxM=reverse(t,Xh,F,Kmid+0*K);
XxH=reverse(t,Xh,F,Khigh+0*K);

E=XxH(:,1)-XxL(:,1);

span=2;
epslocal=zeros(length(t),1);
for k=1:length(t)
    inds=max(1,k-span):length(t);
    %Without zero mean, you lose the nuts for the trees
    %epslocal(k)=(E(inds)-mean(E(inds)))\(XxM(inds,1)-mean(XxM(inds,1)));
    temp=[E(inds) ones(length(inds),1)]\XxM(inds,1);
    epslocal(k)=temp(1);
end
%epslocal=(epslocal(:,1)*Khigh+epslocal(:,2)*Klow)./sum(epslocal,2)

realdot=dot(XxM(:,1),Xd(:,1)/norm(Xd(:,1)))
epsdot=dot(XxL(:,1),E)
epsglobal=E\XxM(:,1)

epsComp=(XxL(:,1)-mean(XxL(:,1)))\(E-mean(E))


figure(3)
clf
subplot(2,1,1)
hold on

plot(t,Xd(:,1),'g')
plot(t,Xh(:,1),'k')
plot(t,Ksmooth/10,'m')

plot(t,XxL(:,1),'ro','markersize',2)
plot(t,XxM(:,1),'rv','markersize',2)
plot(t,XxH(:,1),'r^','markersize',2)
plot(t,XxL(:,1)+epsComp*E,'color',[.8 .5 .2])
plot(t,XxL(:,1)+epslocal.*E,'color',[0 .2 .6])

title('You can recover the twue intent by overestimation correction')
ylabel('Displacement in X')
xlim([0 1.5])
legend('Ground Truth Intention','Hand (Stiffness=magenta)','Time-varying Stiffness/10',['Extraction with estimated stiffness (Stiffness=',num2str(Klow),')'],['Extraction with estimated stiffness (Stiffness=',num2str(Kmid),')'],['Extraction with estimated stiffness (Stiffness=',num2str(Khigh),')'],'Constant Corrected Stiffness','Locally Corrected Stiffness')

subplot(2,1,2)
hold on

%plot(t,epsglobal+0*t,'color',[.8 .5 .2])
%plot(t,epslocal,'color',[0 .2 .6])
plot(t,Xd(:,1)-Xh(:,1),'g')
plot(t,10*(XxH(:,1)-XxM(:,1)),'r')

%plot(t,epslocal(:,1),'b')
%plot(t,(XxM(:,1)-Xd(:,1)),'color',[.7 .4 .4])
%legend('Difference of Estimates','Ground Truth Intention','Constant-Corrected Extraction','Time-varying Stiffness Estimate Extraction')
title(['Solved Stiffness Overestimation: Eps = ',num2str(epsglobal*(Khigh-Klow))])
xlabel('Time, Seconds')
ylabel('Computed scale')

figure(4)
clf
hold on
%plot(t,epsdot*E,'k')
%plot(t,(XxM(:,1)-Xd(:,1)),'color',[.7 .4 .4])
plot(t,0*t,'g')
plot(t,Xh(:,1)-Xd(:,1),'k')
plot(t,XxL(:,1)-Xd(:,1),'ro','markersize',2)
plot(t,XxM(:,1)-Xd(:,1),'rv','markersize',2)
plot(t,XxH(:,1)-Xd(:,1),'r^','markersize',2)
plot(t,XxH(:,1)-XxL(:,1),'c')
plot(t,XxL(:,1)+epsComp*E-Xd(:,1),'r')
plot(t,XxM(:,1)-epslocal.*E-Xd(:,1),'color',[0 .2 .6])

eh=XxL(:,1)-Xh(:,1);
he=E\eh

%plot(t,XxL(:,1)+he*E-Xd(:,1),'color',[.5 .1 .8])


legend({'Zero Error';
    'Hand Error';
    ['Extraction Error (Stiffness=',num2str(Klow),')'];
    ['Extraction Error (Stiffness=',num2str(Kmid),')'];
    ['Extraction Error (Stiffness=',num2str(Khigh),')'];
    'Unscaled Error Estimate';
    'Correction with Linear Regression';
    ['Correction with Local Linear Regression, Span = ',num2str(span)]})
xlabel('Time, sec')
ylabel('Error, X')

