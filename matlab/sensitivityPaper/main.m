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


%% Step 3a: Set up data for mass simulation - Burdet
clear OAT
l1nom=.33;
l2nom=.34;
weightnom=175; %lbs
massnom=weightnom*.4535; %kg
shouldernom=[0 .45];

%Winters (1990)
lc1nom=.436*l1nom;
lc2nom=.682*l2nom;

m1nom=.028*massnom;
m2nom=.022*massnom;

%Shoulder location
x0nom=origin+shouldernom; %Shoulder is measured in room coordinates relative to the workspace center

kp0nom=[10.8 2.83; 2.51 8.67];
kp1nom=[3.18 2.15; 2.34 6.18];
kpkdrationom=1/12;
reflexrationom=1/50;
kpkdreflexrationom=2;

checkratio=1+[-.3 -.2 -.1 -.05 0 .05 .1 .2 .3];
checkruler=[-3 -1 -1 -.5 0 .5 1 2 3]/100;

varyme=ones(18*length(checkratio)+4*length(checkruler),22);
svm1=size(varyme,1);
variedcol=zeros(svm1,1);
ind=0;
for k=[1 2 21 22]
    for kk=1:length(checkruler)
        ind=ind+1;
        varyme(ind,:)=[0 0 ones(1,18) 0 0];
        varyme(ind,k)=checkruler(kk);
        variedcol(ind)=k;
    end
end

for k=3:20
    for kk=1:length(checkratio)
        ind=ind+1;
        varyme(ind,:)=[0 0 ones(1,18) 0 0];
        varyme(ind,k)=checkratio(kk);
        variedcol(ind)=k;
    end
end

tic
figure(4)
clf
hold on
for N=1:svm1
    v=varyme(N,:);
    l1=l1nom+v(1);
    l2=l2nom+v(2);
    lc1=lc1nom*v(3); %Fair because multiplication is commutative
    lc2=lc2nom*v(4);
    m1=m1nom*v(5);
    m2=m2nom*v(6);
    I1=m1*(.322*v(7)*l1)^2;
    I2=m2*(.468*v(8)*l2)^2;
    kp0=v(17)*kp0nom.*[v(9) v(10); v(11) v(12)];
    kp1=v(17)*kp1nom.*[v(13) v(14); v(15) v(16)];
    kpkdratio=kpkdrationom*v(18);
    reflexratio=reflexrationom*v(19);
    kpkdreflexratio=kpkdreflexrationom*v(20);
    x0=x0nom+[v(21) v(22)];
    OAT(N).y=extract(t,xvaf,'reflex');
    plot(OAT(N).y(:,1),OAT(N).y(:,2))
    axis equal
    drawnow
    [N/svm1 toc/N ((svm1/N-1)*(toc))/60]
end

peakV=sqrt(max(vecmag2(xvaf(:,3:4))));
for k=1:length(OAT)
    OAT(k).mue=mean(vecmag(yex(:,1:2)-OAT(k).y(:,1:2)))*1000;
    OAT(k).mueV=mean(vecmag(yex(:,3:4)-OAT(k).y(:,3:4)))/peakV;
end


upper=max([OAT.mue]);
upperV=max([OAT.mueV]);

names={'l1','l2','lc1','lc2','m1','m2','I1','I2','kp0_{11}','kp0_{12}','kp0_{21}','kp0_{22}','kp1_{11}','kp1_{12}','kp1_{21}','kp1_{22}','kpgain','kpkd ratio','reflex ratio','kpkd reflex ratio','P0_x','P0_y'};

figure(667)
clf

figure(668)
clf

for k=1:22
    inds=find(variedcol==k);
    vars=varyme(inds,k);
    vals=[OAT(inds).mue];
    for i=inds'
        OAT(i).variation=varyme(i,k);
        OAT(i).name=names{k};
        OAT(i).col=k;
    end
    figure(667)
    subplot(8,3,k)
    plot(vars,vals,'k')
    ylabel(names{k})
    ylim([0 upper])
    figure(668)
    subplot(8,3,k)
    plot(vars,[OAT(inds).mueV],'k')
    ylabel(names{k})
    ylim([0 upperV])
end

figure(667)
suplabel('Burdet','t')
suplabel('MUE in Position, mm','y')
figure(668)
suplabel('Burdet','t')
suplabel('MUE in Velocity, fraction peak velocity','y')

%% Step 3b: Set up data for mass simulation - Shad & Muss
clear OATSM
l1nom=.33;
l2nom=.34;
weightnom=175; %lbs
massnom=weightnom*.4535; %kg
shouldernom=[0 .45];

%Winters (1990)
lc1nom=.436*l1nom;
lc2nom=.682*l2nom;

m1nom=.028*massnom;
m2nom=.022*massnom;

%Shoulder location
x0nom=origin+shouldernom; %Shoulder is measured in room coordinates relative to the workspace center

kpgain=1;
kpnom=[15 6;6 16];
kdnom=[2.3 .09; .09 2.4];

checkratio=1+[-.3 -.2 -.1 -.05 0 .05 .1 .2 .3];
checkruler=[-3 -1 -1 -.5 0 .5 1 2 3]/100;

varyme=ones(15*length(checkratio)+4*length(checkruler),19);
svm1=size(varyme,1);
variedcol=zeros(svm1,1);
ind=0;
for k=[1 2 18 19]
    for kk=1:length(checkruler)
        ind=ind+1;
        varyme(ind,:)=[0 0 ones(1,15) 0 0];
        varyme(ind,k)=checkruler(kk);
        variedcol(ind)=k;
    end
end

for k=3:17
    for kk=1:length(checkratio)
        ind=ind+1;
        varyme(ind,:)=[0 0 ones(1,15) 0 0];
        varyme(ind,k)=checkratio(kk);
        variedcol(ind)=k;
    end
end

tic
figure(3)
clf
hold on
for N=1:svm1
    v=varyme(N,:);
    l1=l1nom+v(1);
    l2=l2nom+v(2);
    lc1=lc1nom*v(3); %Fair because multiplication is commutative
    lc2=lc2nom*v(4);
    m1=m1nom*v(5);
    m2=m2nom*v(6);
    I1=m1*(.322*v(7)*l1)^2;
    I2=m2*(.468*v(8)*l2)^2;
    kp=v(17)*kpnom.*[v(9) v(10); v(11) v(12)];
    kd=v(17)*kdnom.*[v(13) v(14); v(15) v(16)];
    x0=x0nom+[v(18) v(19)];
    OATSM(N).y=extract(t,xvafsm,@armdynamics_inverted);
    plot(OATSM(N).y(:,1),OATSM(N).y(:,2))
    axis equal
    drawnow
    [N/svm1 toc/N ((svm1/N-1)*(toc))/60]
end

peakV=sqrt(max(vecmag2(xvafsm(:,3:4))));

for k=1:length(OATSM)
    OATSM(k).mue=mean(vecmag(yexsm(:,1:2)-OATSM(k).y(:,1:2)))*1000;
    OATSM(k).mueV=mean(vecmag(yexsm(:,3:4)-OATSM(k).y(:,3:4)))/peakV;
end


upper=max([OATSM.mue]);
upperV=max([OATSM.mueV]);

names={'l1','l2','lc1','lc2','m1','m2','I1','I2','kp11','kp12','kp21','kp22','kd11','kd12','kd21','kd22','kpgain','P0_x','P0_y'};

figure(665)
clf
figure(666)
clf
for k=1:19
    inds=find(variedcol==k);
    vars=varyme(inds,k);
    vals=[OATSM(inds).mue];
    for i=inds'
        OATSM(i).variation=varyme(i,k);
        OATSM(i).name=names{k};
        OATSM(i).col=k;
    end
    figure(665)
    subplot(7,3,k)
    plot(vars,vals,'k')
    ylabel(names{k})
    ylim([0 upper])
    figure(666)
    subplot(7,3,k)
    plot(vars,[OATSM(inds).mueV],'k')
    ylabel(names{k})
    ylim([0 upperV])
end

figure(665)
suplabel('ShadMuss','t')
suplabel('MUE in Position, mm','y')
figure(666)
suplabel('ShadMuss','t')
suplabel('MUE in Velocity, fraction peak velocity','y')

save('OATs.mat','OAT','OATSM')

%% Step 4a: Monte Carlo Variance estimation - Burdet

% First step, set up nominal values.
if exist('simSobol.mat','file')
    return
end

l1nom=.33;
l2nom=.34;
weightnom=175; %lbs
massnom=weightnom*.4535; %kg
shouldernom=[0 .45];

%Winters (1990)
lc1nom=.436*l1nom;
lc2nom=.682*l2nom;

m1nom=.028*massnom;
m2nom=.022*massnom;

%rog of gyration numbers from winters, rog=sqrt(I/m)
I1nom=m1nom*(.322*l1nom)^2;
I2nom=m2nom*(.468*l2nom)^2;

%Shoulder location
x0=origin+shouldernom; %Shoulder is measured in room coordinates relative to the workspace center

kpgain=1;
kp0nom=[10.8 2.83; 2.51 8.67];
kp1nom=[3.18 2.15; 2.34 6.18];

% Next, which parameters are really independent? Certainly not the
% terms inside a stiffness matrix...
% l1 l2 lc1 lc2 m1 m2 i1 i2 kp0 kp1... but really you want to fuzz up
% winters' numbers.

p=sobolset(20,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
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

save('sobolSM.mat','saltelliA','saltelliB','saltelliAB');

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

varY=var([saltelliA(:).mue]);
varYv=var([saltelliA(:).mueV]);

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