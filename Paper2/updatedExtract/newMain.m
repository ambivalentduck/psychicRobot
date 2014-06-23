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
K1nom=15;
K2nom=16;
Koffdiagnom=6;

% K1nom=8;
% K2nom=11;
% Koffdiagnom=3;

K1diag=K1nom;
K2diag=K2nom;
Koffdiag=Koffdiagnom;

measuredVals=qt;

kN=8;
minKr=.5;
maxKr=1.5;

mFfilt=zeros(length(t),4);

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

R2k1=t;
m1d=t;
for k=1:length(t)
    temp=[klist 0*klist+1]\(linme(k,:)');
    m1d(k)=temp(1);
    R2k1(k)=corr(linme(k,:)',klist).^2;
end
mFfilt(:,1)=m1d;
subplot(6,2,2)
plot(t,R2k1,'b.')
ylabel('R^2')
title('Quality and Values of Regression of K estimate onto Extraction')
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
mFfilt(:,2)=m2d;
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
mFfilt(:,3)=mod1;
mFfilt(:,4)=mod2;
subplot(6,2,10)
plot(t,R21,'b.',t,R22,'k.')
ylabel('R^2')
subplot(6,2,12)
plot(t,mod1,'b.',t,mod2,'k.')
ylabel('Time-varying Slope')
xlabel('Time, s')

%% We have some time-varying m times Ffilt. Now for K

figure(2)
clf
K1diag=K1nom;
K2diag=K2nom;
Koffdiag=Koffdiagnom;

inds=find((t<.6)&(t>0));
[T,Q]=ode45(@armdynamics_inverted,t,qt(1,1:4));

%K1
subplot(6,2,[1 3])
hold on
x=q2x(Q);
plot(x(:,1),x(:,2),'r')
x=q2x([Q(:,1)-K1nom*mFfilt(:,1) Q(:,2)]);
plot(x(:,1),x(:,2),'c')
plot(xvaf(:,1),xvaf(:,2),'k')

ycalc=Q(inds,1)+K1nom*mFfilt(inds,1);
mF=mFfilt(inds,1);
mF=mF-mean(mF);
K1=dot(ycalc-mean(ycalc),mF)/dot(mF,mF)
x=q2x([Q(:,1)+(K1nom-K1)*mFfilt(:,1) Q(:,2)]);
plot(x(:,1),x(:,2),'b')
axis equal

subplot(6,2,2)
hold on
plot(t,qt(:,1),'k')
plot(t,Q(:,1),'r')
plot(t,mFfilt(:,1)*10,'b')
title(['K1 calc = ',num2str(K1)])

%K2
subplot(6,2,[5 7])
hold on
plot(xvaf(:,1),xvaf(:,2),'k')
x=q2x(Q);
plot(x(:,1),x(:,2),'r')
x=q2x([Q(:,1) Q(:,2)-K2nom*mFfilt(:,2)]);
plot(x(:,1),x(:,2),'c')

ycalc=Q(inds,2)+K2nom*mFfilt(inds,2);
mF=mFfilt(inds,2);
mF=mF-mean(mF);
K2=dot(ycalc-mean(ycalc),mF)/dot(mF,mF)
x=q2x([Q(:,1) Q(:,2)+(K2nom-K2)*mFfilt(:,2)]);
plot(x(:,1),x(:,2),'b')
axis equal

subplot(6,2,6)
hold on
plot(t,qt(:,1),'k')
plot(t,Q(:,1),'r')
plot(t,mFfilt(:,2)*10,'b')
title(['K1 calc = ',num2str(K2)])

%KOD
subplot(6,2,[9 11])
hold on
plot(xvaf(:,1),xvaf(:,2),'k')
x=q2x(Q);
plot(x(:,1),x(:,2),'r')
x=q2x([Q(:,1)-Koffdiagnom*mFfilt(:,3) Q(:,2)-Koffdiagnom*mFfilt(:,4)]);
plot(x(:,1),x(:,2),'c')

ycalc=[Q(inds,1)+Koffdiagnom*mFfilt(:,3); Q(inds,2)+Koffdiagnom*mFfilt(:,4)];
mF=[mFfilt(inds,3); mFfilt(inds,4)];
mF=mF-mean(mF);
KOD=dot(ycalc-mean(ycalc),mF)/dot(mF,mF)
x=q2x([Q(:,1) Q(:,2)+(K2nom-K2)*mFfilt(:,2)]);
plot(x(:,1),x(:,2),'b')
axis equal

subplot(6,2,6)
hold on
plot(t,qt(:,1),'k')
plot(t,Q(:,1),'r')
plot(t,mFfilt(:,2)*10,'b')
title(['K1 calc = ',num2str(K2)])

%% Use all of that craziness


K1diag=K1;
K2diag=K2;
Koffdiag=KOD;
[T,Q]=ode45(@armdynamics_inverted,t(inds),qt(1,1:4));
figure(3)
clf
hold on
x=q2x(Q);
plot(xvaf(:,1),xvaf(:,2),'k')
plot(x(:,1),x(:,2),'r')
axis equal


%% And we're done
cleanup