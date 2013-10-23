clc
clear all

global measuredVals measuredTime l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kp kd kpgain kpkdratio reflexratio kpkdreflexratio

%% Step 0: Set physical parameters

LEFT=.281;
RIGHT=-.3951;
TOP=.18;
BOTTOM=.68;
origin=[(LEFT+RIGHT)/2,(TOP+BOTTOM)/2]; %Might as well duplicate the physical workspace...why not?

%Use my own physical parameters because...why not? Also, no IRB/justification
l1=.33;
l2=.34;
weight=175; %lbs
mass=weight*.4535; %kg
shoulder=[0 .45];

%Winters (1990)
lc1=.436*l1;
lc2=.682*l2;

m1=.028*mass;
m2=.022*mass;

%rog of gyration numbers from winters, rog=sqrt(I/m)
I1=m1*(.322*l1)^2;
I2=m2*(.468*l2)^2;

%Shoulder location
x0=origin+shoulder; %Shoulder is measured in room coordinates relative to the workspace center

kp0=[10.8 2.83; 2.51 8.67];
kp1=[3.18 2.15; 2.34 6.18];
kp=[15 6;6 16];
kd=[2.3 .09; .09 2.4];
kpkdratio=1/12;
reflexratio=1/50;
kpkdreflexratio=2;


%% Step 1: Set up a basic kicked movement and forward simulate.

t=0:.005:2;
coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
[x,v,a]=minjerk(coeff,t);
x=x';
v=v';
a=a';

i=find(t>=.7);
x(i,1)=x(i(1)-1,1);
x(i,2)=x(i(1)-1,2);
v(i,:)=0;
a(i,:)=0;

f=zeros(length(t),2);

i=find((t>=.1)&(t<=.15));
f(i,2)=15;

xvaf=[x v a f];

measuredVals=xvaf;
measuredTime=t;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q,l1,l2);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)',l1,l2);
    torque=-fJq'*xvaf(k,7:8)';
    measuredVals(k,:)=[q' qdot' qddot' torque'];
end

q0=measuredVals(1,1:4);

%forward simulate
kpgain=1;
[T,X]=forwardReflexHelper(t,q0);
[T,XSM]=ode45(@armdynamicsShadMuss,t,q0);

y=zeros(length(T),6);
ysm=y;

for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2),l1,l2)*X(k,3:4)')';
    ysm(k,1:2)=fkin(XSM(k,1:2));
    ysm(k,3:4)=(fJ(XSM(k,1:2),l1,l2)*XSM(k,3:4)')';
end

gT=gradient(t)';
y(:,5)=gradient(y(:,3))./gT;
y(:,6)=gradient(y(:,4))./gT;
ysm(:,5)=gradient(ysm(:,3))./gT;
ysm(:,6)=gradient(ysm(:,4))./gT;


figure(1)
clf
hold on
plot(x(:,1),x(:,2),'b',y(:,1),y(:,2),'k.')
quiver(y(:,1),y(:,2),f(:,1),f(:,2),'b')
plot(ysm(:,1),ysm(:,2),'k^')
axis equal

%% Step 2: Extract to demonstrate accuracy

yex=extract(t,[y f],'reflex');
yexsm=extract(t,[ysm f],@armdynamics_inverted);
plot(yex(:,1),yex(:,2),'r.')
plot(yexsm(:,1),yexsm(:,2),'r^')
ycross=extract(t,[y f],@armdynamics_inverted);
plot(ycross(:,1),ycross(:,2),'r*')
plot(x(:,1),x(:,2),'b')

legend('Intent','Forward Sim','Forces','Shad&Muss','Extracted Intent','Extracted ShadMuss','ShadMuss Extract from Burdet')

xvaf=[y f];
xvafsm=[ysm f];
save('baselines.mat','yex')


%% Step 3a: OAT - Burdet

names=paramsPopulator('names');
dat=paramsPopulator('Burdet');
f=find(dat(:,3));

lf=length(f);
for k=1:lf
    params=paramsPopulator(f(k));
    OAT(k,:)=repeatedSim(params,t,xvaf,'reflexes',0,0);
    for kk=1:size(params,1)
        OAT(k,kk).name=names(f(k));
        OAT(k,kk).val=params(kk,f(k));
    end
end

%% Step 3b: OAT - ShadMuss

dat=paramsPopulator('shadmuss');
f=find(dat(:,3));

lf=length(f);
for k=1:lf
    params=paramsPopulator(f(k));
    OATSM(k,:)=repeatedSim(t,xvaf,'reflexes',0,0);
    for kk=1:size(params,1)
        OATSM(k,kk).name=names(f(k));
        OATSM(k,kk).val=params(kk,f(k));
    end
end

save('OAT_KICK.mat','OAT','OATSM')


%% Step 4a: Monte Carlo Variance estimation - Burdet

p=sobolset(,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

varied=p(1:1000,:); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.
%The number after p controls the number of points. Remember, this N*10 is the number of individual sims done.

consistent=varied(:,1:10);
individual=varied(:,11:20);

%Notice that you're actually slanting things here by giving everything but
%kp0 and kp1 a 10% total range vs kps which have a 20% total range
range=.1*[l1 l2 lc1 lc2 m1 m2 I1 I2]; %Have to vary about nominal values to avoid correlated variance
% Need to fix this.
% converted=norminv(uniform from sobol,mean, sd);
% Probably need to quietly saturate to +/- 3 sd. Histogram of doing that looks fine.
%l1, l2 mean error 0, sd 2 mm
%mass mean error +1.2kg, sd 3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)

sc1=size(consistent,1);

tic
for N=1:size(consistent,1)
    v=consistent(N,:);
    v(1:8)=(v(1:8)-.5).*range;
    v(9:10)=(v(9:10)-.5)*.4+1;
    l1=norminv(v(1),l1nom,.01); %Measurement error on order of mm
    l2=norminv(v(2),l2nom,.01);
    lc1=norminv(v(3),.436,.01/(.436*l1nom))*l1; %choose sd ~1cm, Dempster does NOT have this
    lc2=norminv(v(4),.682,.01/(.682*l2nom))*l2; %choose sd ~1cm
    mass=norminv(v(5),massnom,3.1); %mass self-report error SD=3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)
    m1=norminv(v(6),.028,.0029)*mass; % Dempster 1955
    m2=norminv(v(7),.022,.00248)*mass; % Dempster 1955
    I1=m1*(norminv(v(8),.322,.322/20)*l1)^2; % ~5% SD
    I2=m2*(norminv(v(9),.468,.468/20)*l2)^2; % ~5% SD
    kp=kpnom*norminv(v(10),1,.15); % SD is actually quite large
    kd=kdnom*norminv(v(11),1,.15); % SD is actually quite large
    simmedA(N).y=extract(t,xvaf,'reflex');

    v=individual(N,:);
    v(1:8)=(v(1:8)-.5).*range;
    v(9:10)=(v(9:10)-.5)*.4+1;
    l1=l1nom+v(1);
    l2=l2nom+v(2);
    lc1=lc1nom+v(3);
    lc2=lc2nom+v(4);
    m1=m1nom+v(5);
    m2=m2nom+v(6);
    I1=I1nom+v(7);
    I2=I2nom+v(8);
    kp0=kp0nom*v(9);
    kp1=kp1nom*v(10);
    simmedB(N).y=extract(t,xvaf,'reflex');

    for col=1:10
        v=consistent(N,:);
        v(col)=individual(N,col); %Overwrite just a single value at a time.
        v(1:8)=(v(1:8)-.5).*range;
        v(9:10)=(v(9:10)-.5)*.4+1;
        l1=l1nom+v(1);
        l2=l2nom+v(2);
        lc1=lc1nom+v(3);
        lc2=lc2nom+v(4);
        m1=m1nom+v(5);
        m2=m2nom+v(6);
        I1=I1nom+v(7);
        I2=I2nom+v(8);
        kp0=kp0nom*v(9);
        kp1=kp1nom*v(10);
        simmedAB(N,col).y=extract(t,xvaf,'reflex');
    end
    [N/sc1 toc/N ((sc1/N-1)*(toc))/60]
end

%% Plot the results of the quasi-Monte Carlo analysis

load('simSobol.mat')

names={'l1','l2','lc1','lc2','m1','m2','I1','I2','kp0','kp1'};

mueA=zeros(length(simmedA),1);
mueB=mueA;
mueAB=zeros(length(simmedA),10);

for k=1:length(simmedA)
    mueA(k)=mean(vecmag(yex(:,1:2)-simmedA(k).y(:,1:2)))*1000;
    mueB(k)=mean(vecmag(yex(:,1:2)-simmedB(k).y(:,1:2)))*1000;
    for col=1:10
        mueAB(k,col)=mean(vecmag(yex(:,1:2)-simmedAB(k,col).y(:,1:2)))*1000;
    end
end

EnxVx=zeros(10,1);
VxEnx=EnxVx;

mueABm=mueAB;
for col=1:10
    EnxVx(col)=1/length(simmedA)*sum(mueB.*(mueAB(:,col)-mueA));
    VxEnx(col)=1/(2*length(simmedA))*sum((mueA-mueAB(:,col)).^2);
end

%varA=var(mueA);
varA=var(mueAB(:))

S=VxEnx/varA
STi=EnxVx/varA

figure(237)
clf
hold on
plot([1-.15 10.15],[0 0],'m')
for col=1:10
    plot(col*ones(length(simmedA),1)-.15,mueA-mueAB(:,col),'k.',col-.15,VxEnx(col),'rx')
    plot(col*ones(length(simmedA),1)-.05,mueB.*(mueAB(:,col)-mueA),'b.',col-.05,EnxVx(col),'rx')
end
ylabel('Mean Unsigned Error, millimeters')
set(gca,'xtick',1:10)
set(gca,'xticklabel',names)

%% Monte Carlo Variance estimation - ShadMuss

if exist('sobolSM.mat','file')
    return
end

l1nom=.33;
l2nom=.34;
weightnom=175; %lbs
massnom=weightnom*.4535; %kg
shouldernom=[0 .45];

%Shoulder location
x0=origin+shouldernom; %Shoulder is measured in room coordinates relative to the workspace center

kpgain=1;
kpnom=[15 6;6 16];
kdnom=[2.3 .09; .09 2.4];

p=sobolset(22,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

varied=p(1:1000,:); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.

%The number after p controls the number of points. Remember, this N*10 is the number of individual sims done.

A=varied(:,1:11);
B=varied(:,12:22);

tic
for N=1:svm1
    v=A(N,:);
    v(v<.001)=.001; %Avoid having to deal with infinity and other impossible values
    v(v>.999)=.999;
    l1=norminv(v(1),l1nom,.001); %Measurement error on order of mm
    l2=norminv(v(2),l2nom,.001);
    lc1=norminv(v(3),.436,.01/(.436*l1nom))*l1; %choose sd ~1cm
    lc2=norminv(v(4),.682,.01/(.682*l2nom))*l2; %choose sd ~1cm
    mass=norminv(v(5),massnom,3.1); %mass self-report error SD=3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)
    m1=norminv(v(6),.028,.028/20)*mass; % ~5% SD
    m2=norminv(v(7),.022,.022/20)*mass; % ~5% SD
    I1=m1*(norminv(v(8),.322,.322/20)*l1)^2; % ~5% SD
    I2=m2*(norminv(v(9),.468,.468/20)*l2)^2; % ~5% SD
    kp=kpnom*norminv(v(10),1,.15); % SD is actually quite large
    kd=kdnom*norminv(v(11),1,.15); % SD is actually quite large
    saltelliA(N).y=extract(t,xvafsm,@armdynamics_inverted);

    v=B(N,:);
    v(v<.001)=.001; %Avoid having to deal with infinity and other impossible values
    v(v>.999)=.999;
    l1=norminv(v(1),l1nom,.001); %Measurement error on order of mm
    l2=norminv(v(2),l2nom,.001);
    lc1=norminv(v(3),.436,.01/(.436*l1nom))*l1; %choose sd ~1cm
    lc2=norminv(v(4),.682,.01/(.682*l2nom))*l2; %choose sd ~1cm
    mass=norminv(v(5),massnom,3.1); %mass self-report error SD=3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)
    m1=norminv(v(6),.028,.028/20)*mass; % ~5% SD
    m2=norminv(v(7),.022,.022/20)*mass; % ~5% SD
    I1=m1*(norminv(v(8),.322,.322/20)*l1)^2; % ~5% SD
    I2=m2*(norminv(v(9),.468,.468/20)*l2)^2; % ~5% SD
    kp=kpnom*norminv(v(10),1,.15); % SD is actually quite large
    kd=kdnom*norminv(v(11),1,.15); % SD is actually quite large
    saltelliB(N).y=extract(t,xvafsm,@armdynamics_inverted);

    for k=1:11
        v=A(N,:);
        v(k)=B(N,k);
        v(v<.001)=.001; %Saturate to avoid having to deal with infinity and other impossible values
        v(v>.999)=.999;

        l1=norminv(v(1),l1nom,.001); %Measurement error on order of mm
        l2=norminv(v(2),l2nom,.001);
        lc1=norminv(v(3),.436,.01/(.436*l1nom))*l1; %choose sd ~1cm
        lc2=norminv(v(4),.682,.01/(.682*l2nom))*l2; %choose sd ~1cm
        mass=norminv(v(5),massnom,3.1); %mass self-report error SD=3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)
        m1=norminv(v(6),.028,.028/20)*mass; % ~5% SD
        m2=norminv(v(7),.022,.022/20)*mass; % ~5% SD
        I1=m1*(norminv(v(8),.322,.322/20)*l1)^2; % ~5% SD
        I2=m2*(norminv(v(9),.468,.468/20)*l2)^2; % ~5% SD
        kp=kpnom*norminv(v(10),1,.15); % SD is actually quite large
        kd=kdnom*norminv(v(11),1,.15); % SD is actually quite large
        saltelliAB(N,k).y=extract(t,xvafsm,@armdynamics_inverted);
    end

    [N/svm1 toc/N ((svm1/N-1)*(toc))/60]
end

peakV=sqrt(max(vecmag2(xvafsm(:,3:4))));

for k=1:length(saltelliA)
    saltelliA(k).mue=mean(vecmag(yexsm(:,1:2)-saltelliA(k).y(:,1:2)))*1000;
    saltelliA(k).mueV=mean(vecmag(yexsm(:,3:4)-saltelliA(k).y(:,3:4)))/peakV;

    saltelliB(k).mue=mean(vecmag(yexsm(:,1:2)-saltelliB(k).y(:,1:2)))*1000;
    saltelliB(k).mueV=mean(vecmag(yexsm(:,3:4)-saltelliB(k).y(:,3:4)))/peakV;
    for kk=1:11
        saltelliAB(k,kk).mue=mean(vecmag(yexsm(:,1:2)-saltelliAB(k,kk).y(:,1:2)))*1000;
        saltelliAB(k,kk).mueV=mean(vecmag(yexsm(:,3:4)-saltelliAB(k,kk).y(:,3:4)))/peakV;
    end
end

save('sobolSM.mat','saltelliA','saltelliB','saltelliAB','yexsm');

%% Plot outcome of sensitivity analysis on ShadMuss

load sobolSM.mat

EnxVx=zeros(11,1);
VxEnx=EnxVx;
EnxVxV=zeros(11,1);
VxEnxV=EnxVx;

fA=[saltelliA.mue];
fB=[saltelliB.mue];
fvA=[saltelliA.mueV];
fvB=[saltelliB.mueV];

varY=var([saltelliAB(:).mue]);
varYv=var([saltelliAB(:).mueV]);

figure(237)
clf
hold on
figure(238)
clf
hold on
for col=1:11
    fAB=[saltelliAB(:,col).mue];
    fvAB=[saltelliAB(:,col).mueV];
    EnxVx(col)=1/length(fA)*sum(fB.*(fAB-fA));
    VxEnx(col)=1/(2*length(fA))*sum((fA-fAB).^2);
    EnxVxV(col)=1/length(fA)*sum(fvB.*(fvAB-fvA));
    VxEnxV(col)=1/(2*length(fA))*sum((fvA-fvAB).^2);

    figure(237)
    plot(col*ones(length(fA),1)-.1,(fA-fAB)/varY,'k.',col-.1,VxEnx(col)/varY,'rx')
    plot(col*ones(length(fA),1)+.1,fB.*(fAB-fA)/varY,'b.',col+.1,EnxVx(col)/varY,'rx')
    figure(238)
    %plot(col*ones(length(fA),1)-.1,(fvA-fvAB)/varYv,'k.',col-.1,VxEnxV(col)/varYv,'rx')
    %plot(col*ones(length(fA),1)+.1,fvB.*(fvAB-fvA)/varYv,'b.',col+.1,EnxVxV(col)/varYv,'rx')
    plot(col-.1,VxEnxV(col)/varYv,'kx')
    plot(col+.1,EnxVxV(col)/varYv,'bx')
end
figure(237)
plot([1-.15 11.15],[0 0],'Color',[0 .7 0])
ylabel('Mean Unsigned Error, fraction total')
set(gca,'xtick',1:11)
names={'l1','l2','lc1','lc2','mass','m1','m2','I1','I2','kp','kd'};
set(gca,'xticklabel',names)
title('Position')
ylim([-2 2])

figure(238)
plot([1-.1 11.1],[0 0],'Color',[0 .7 0])
ylabel('Mean Unsigned Error, fraction total')
set(gca,'xtick',1:11)
names={'l1','l2','lc1','lc2','mass','m1','m2','I1','I2','kp','kd'};
set(gca,'xticklabel',names)
title('Velocity')

S=VxEnx/varY
STi=EnxVx/varY

SV=VxEnxV/varYv
STiV=EnxVxV/varYv

figure(239)
clf
bar(1:length(S),[S SV STi STiV])
names={'l1','l2','lc1','lc2','mass','m1','m2','I1','I2','kp','kd'};
set(gca,'xticklabel',names)
ylabel('Mean Unsigned Error, fraction total')
legend('S Pos','S Vel','ST Pos','ST Vel')