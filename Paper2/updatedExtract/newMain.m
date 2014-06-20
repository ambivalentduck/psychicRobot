%% Housekeeping

clc
clear all

load('pulse1.mat');
global measuredVals measuredTime K1diag K2diag Koffdiag

set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass);

%% Snag an early pulse and convert xvaf to q and torque

dcats=[trials.disturbcat];
f=find(dcats==1);
N=f(4);

xvaf=[trials(N).x trials(N).v trials(N).a trials(N).f];

onset=find(vecmag(trials(N).v)>.05,1,'first');
start=max(onset-35,1);


qt=xvaf2arm(xvaf(start:end,:));
t=trials(N).t(start:end);
t=t-t(1);
measuredTime=t;

%% Verify that we're still time-locally linear in each of the K parameters

%Params from Shad&Muss
K1diag=15;
K1nom=15;
clf

K2diag=16;
K2nom=16;
Koffdiag=6;
Koffdiagnom=6;

measuredVals=qt;

kN=8;
minKr=.5;
maxKr=1.5;

figure(1)
clf
subplot(6,2,[1 3])
hold on
plot(xvaf(:,1),xvaf(:,2),'k.')

linme=zeros(length(t),kN);
klist=zeros(kN,1);
for k=1:kN
    K1diag=K1nom*(minKr+(maxKr-minKr)*(k-1)/(kN-1));
    klist(k)=K1diag;
    [T,Q]=ode45(@armdynamics_inverted,t,qt(1,1:4));
    linme(:,k)=Q(:,1);
    x=q2x(Q);
    plot(x(:,1),x(:,2),'r.')
end
K1diag=K1nom;
title('K1 diagonal element is varied.')
axis equal
axis off

R2=t;
m1d=t;
for k=1:length(t)
    temp=[klist 0*klist+1]\(linme(k,:)');
    m1d(k)=temp(1);
    R2(k)=corr(linme(k,:)',klist).^2;
end
subplot(6,2,2)
plot(t,R2,'b.')
ylabel('R^2')
subplot(6,2,4)
plot(t,m1d,'b.')
ylabel('Time-varying Slope')
xlabel('Time, s')

subplot(6,2,[5 7])
hold on
plot(xvaf(:,1),xvaf(:,2),'k.')

linme=zeros(length(t),kN);
klist=zeros(kN,1);
for k=1:kN
    K2diag=K2nom*(minKr+(maxKr-minKr)*(k-1)/(kN-1));
    klist(k)=K2diag;
    [T,Q]=ode45(@armdynamics_inverted,t,qt(1,1:4));
    linme(:,k)=Q(:,2);
    x=q2x(Q);
    plot(x(:,1),x(:,2),'r.')
end
K2diag=K2nom;
title('K2 diagonal element is varied.')
axis equal
axis off

R2=t;
m2d=t;
for k=1:length(t)
    temp=[klist 0*klist+1]\(linme(k,:)');
    m2d(k)=temp(1);
    R2(k)=corr(linme(k,:)',klist).^2;
end
subplot(6,2,6)
plot(t,R2,'k.')
ylabel('R^2')
subplot(6,2,8)
plot(t,m1d,'k.')
ylabel('Time-varying Slope')
xlabel('Time, s')

subplot(6,2,[9 11])
hold on
plot(xvaf(:,1),xvaf(:,2),'k.')

linme1=zeros(length(t),kN);
linme2=zeros(length(t),kN);
klist=zeros(kN,1);
for k=1:kN
    K1diag=K1nom*(minKr+(maxKr-minKr)*(k-1)/(kN-1));
    klist(k)=K1diag;
    [T,Q]=ode45(@armdynamics_inverted,t,qt(1,1:4));
    linme1(:,k)=Q(:,1);
    linme2(:,k)=Q(:,2);
    x=q2x(Q);
    plot(x(:,1),x(:,2),'r.')
end
Koffdiag=Koffdiagnom;
title('Off-diagonal element is varied.')
axis equal
axis off

R21=t;
R22=t;
mod1=t;
mod2=t;
for k=1:length(t)
    temp=[klist 0*klist+1]\(linme1(k,:)');
    mod1(k)=temp(1);
    R21(k)=corr(linme1(k,:)',klist).^2;
    temp=[klist 0*klist+1]\(linme2(k,:)');
    mod2(k)=temp(1);
    R22(k)=corr(linme2(k,:)',klist).^2;
end
subplot(6,2,10)
plot(t,R21,'b.',t,R22,'k.')
ylabel('R^2')
subplot(6,2,12)
plot(t,mod1,'b.',t,mod2,'k.')
ylabel('Time-varying Slope')
xlabel('Time, s')

cleanup