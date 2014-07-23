clc
clear all
close all

global M B

%% Constant K, prove tautology and other housekeeping
tspan=.5;
tinc=.005;
ta=(0:tinc:tspan)/tspan;

xi=[1.1; .1];
xf=[.26; .1];

xd=(xi*ones(size(ta))+(xf-xi)*(10*ta.^3-15*ta.^4+6*ta.^5))';
vd=((xf-xi)*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';
ad=((xf-xi)*(60*ta-180*ta.^2+120*ta.^3)/(tspan^2))';
% tshort=0:.001:tspan;
% xd=([1; 0]*sin(tshort))';
% vd=([1; 0]*cos(tshort))';
% ad=([1; 0]*-sin(tshort))';

t=0:tinc:1.5;
Xd=[[xd vd ad];zeros(length(t)-length(ta),6)];
Xd(length(ta)+1:length(t),1)=Xd(length(ta)+1:length(t),1)+xf(1);
Xd(length(ta)+1:length(t),2)=Xd(length(ta)+1:length(t),2)+xf(2);
[Fb,Fa]=butter(10,3*tinc,'low');
F=[7*filtfilt(Fb,Fa,randn(2*length(t),1)) 7*filter(Fb,Fa,randn(2*length(t),1))];
F=F(end-length(t)+1:end,:);

%F=zeros(length(t),2);
%F((t>.1)&(t<.2),:)=-15;
K=15*ones(length(t),2);

M=.41;
B=2.3;

Xh=forward(t,Xd,F,K);
Xx=reverse(t,Xh,F,K);

figure(1)
clf
subplot(2,1,1)
hold on
plot(Xd(:,1),Xd(:,2),'g')
plot(Xh(:,1),Xh(:,2),'k')
axis equal

subplot(2,1,2)
hold on
plot(t,Xd(:,1),'g')
plot(t,Xh(:,1),'k')
plot(t,Xx(:,1),'r.')
arrow([t' Xh(:,1)],[t' Xh(:,1)]+.05*[0*t' F(:,1)],.6*[1 1 1],.2);

title('Proof of Tautology Being Intact')
xlabel('Time, Seconds')
ylabel('Displacement in X')
legend('Ground Truth Intention','Hand (Stiffness=15)','Intention (Stiffness=15)','Forces')

%% Sobol spanned extractions

sParams=1;
nCandidates=10;
Kml=15.72;
Kspan=.2;

p=sobolset(sParams,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

candidateK=Kml+Kml*Kspan*(p(1:nCandidates,:)); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.

%You know you have some settling time, so figure out what it is.

Xml=reverse(t,Xh,F,Kml+K*0);

Xs=zeros(length(t),nCandidates);
Ys=zeros(length(t),nCandidates);
for k=1:nCandidates
    cand{k}=reverse(t,Xh,F,candidateK(k)+K*0);
    Xs(:,k)=cand{k}(:,1);
    Ys(:,k)=cand{k}(:,2);
    %plot(t,Xs(:,k),'b')
end

eps=zeros(length(t),2);
for k=1:length(t)
    temp=[candidateK 1+0*candidateK]\Xs(k,:)';
    eps(k,1)=temp(1);
    temp=[candidateK 1+0*candidateK]\Ys(k,:)';
    eps(k,2)=temp(1);
end

%% Try to find parameter-error component of intent

ws=zeros(4,2,length(t));

P=eye(4);

%w0=[(15-Kml)*eye(2);eye(2)];
w0=[0*eye(2);eye(2)];

for k=1:length(t)
    %[w,P]=RLS([eps(k,1);eps(k,2); 1],Xd(k,1:2),w,P,50);
    alpha1(k,:)=Xd(k,1:2)-[eps(k,1);eps(k,2); Xml(k,1:2)']'*w0;
    [w,P]=RLS([eps(k,1);eps(k,2); Xml(k,1:2)'],Xd(k,1:2),w,P,1.2);
    alpha2(k,:)=Xd(k,1:2)-[eps(k,1);eps(k,2); Xml(k,1:2)']'*w;
    ws(:,:,k)=w;
end
alpha1-alpha2

figure(2)
clf

subplot(3,1,1)
hold on
plot(Xd(:,1),Xd(:,2),'g','linewidth',2)
plot(Xh(:,1),Xh(:,2),'k','linewidth',2)
plot(Xml(:,1)-(Kml-15)*eps(:,1),Xml(:,2)-(Kml-15)*eps(:,2),'k-')
plot(squeeze(ws(3,1,:)),squeeze(ws(4,2,:)),'r.')
axis equal

subplot(3,1,2)
hold on
plot(t,Xd(:,1),'g','linewidth',2)
plot(t,Xh(:,1),'k','linewidth',2)
plot(t,Xml(:,1)-(Kml-15)*eps(:,1),'k-')
plot(t,squeeze(ws(3,1,:)),'r.')

subplot(3,1,3)
hold on
for k=2:2
    for kk=1:2
        plot(t,squeeze(ws(k,kk,:)),'r')
    end
end
ylim([-5 5])