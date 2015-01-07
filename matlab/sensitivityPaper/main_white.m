clc
clear all

%% Step 1: Set up a basic kicked movement and forward simulate.

t=0:.005:2;
coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
tcalc=t;
tcalc(t>=.7)=.7;
[x,v,a]=minjerk(coeff,tcalc);
x=x';
v=v';
a=a';

blah=find(t>.7);
w=randn(2*length(t),2);
w(length(t)+blah,:)=0; %turn it off when you get close
[Fb,Fa]=butter(10,3*.005,'low');
f=7*filter(Fb,Fa,w);
f=f(length(t)+1:end,:);

xvaf=[x v a f];

measuredVals=xvaf;
measuredTime=t;

[xsim,xsimSM]=forwardSim(paramsPopulator,t,xvaf);

figure(1)
clf
hold on
plot(x(:,1),x(:,2),'b',xsim(:,1),xsim(:,2),'k.')
quiver(xsim(:,1),xsim(:,2),f(:,1),f(:,2),'b')
plot(xsimSM(:,1),xsimSM(:,2),'k^')
axis equal

%% Step 2: Extract to demonstrate accuracy

yex=extract(t,[xsim f],'reflex');
yexsm=extract(t,[xsimSM f],@armdynamics_inverted);
plot(yex(:,1),yex(:,2),'r.')
plot(yexsm(:,1),yexsm(:,2),'r^')
ycross=extract(t,[xsim f],@armdynamics_inverted);
plot(ycross(:,1),ycross(:,2),'m*')
plot(x(:,1),x(:,2),'b')

legend('Intent','Forward Sim','Forces','Shad&Muss','Extracted Intent','Extracted ShadMuss','ShadMuss Extract from Burdet')

xvaf=[xsim f];
xvafsm=[xsimSM f];
drawnow

save('baselines_white.mat','yex','yexsm')

%% Step 3a: OAT - Burdet
if 0
    clear OAT
    
    names=paramsPopulator('names');
    dat=paramsPopulator('burdet');
    f=find(dat(:,3));
    
    params=paramsPopulator(f(1));
    sp1=size(params,1);
    lf=length(f);
    
    vals=zeros(lf,sp1);
    
    for k=1:lf
        params=paramsPopulator(f(k));
        if k==1
            OAT=repeatedSim(params,t,xvaf,'reflex',0,lf);
        else
            OAT(k,:)=repeatedSim(params,t,xvaf,'reflex',0,lf-k+1);
        end
        for kk=1:sp1
            vals(k,kk)=params(kk,f(k));
        end
    end
    
    for k=1:lf
        for kk=1:sp1
            OAT(k,kk).name=names{f(k)};
            OAT(k,kk).val=vals(k,kk);
        end
    end
end

%% Step 3b: OAT - ShadMuss
if 0
    clear OATSM
    dat=paramsPopulator('shadmuss');
    f=find(dat(:,3));
    
    params=paramsPopulator(f(1));
    sp1=size(params,1);
    lf=length(f);
    
    vals=zeros(lf,sp1);
    
    for k=1:lf
        params=paramsPopulator(f(k));
        if k==1
            OATSM=repeatedSim(params,t,xvafsm,@armdynamics_inverted,0,lf);
        else
            OATSM(k,:)=repeatedSim(params,t,xvafsm,@armdynamics_inverted,0,lf-k+1);
        end
        for kk=1:sp1
            vals(k,kk)=params(kk,f(k));
        end
    end
    
    for k=1:lf
        for kk=1:sp1
            OATSM(k,kk).name=names{f(k)};
            OATSM(k,kk).val=vals(k,kk);
        end
    end
    
    save('OAT_KICK_white.mat','OAT','OATSM')
end

%% Step 4a: Set up Monte Carlo Variance estimation

p=sobolset(40,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

varied=p(1:1000,:); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.
%The number after p controls the number of points. Remember, this N*20 is the number of individual sims done.

%Clamp to +/- 3 stdevs to avoid support issues like negative segment masses
varied(varied>.999)=.999;
varied(varied<.001)=.011;

A=varied(:,1:20);
B=varied(:,21:40);

AB=zeros(size(A,1),size(A,2),20);
for k=1:20
    AB(:,:,k)=A;
    AB(:,k,k)=B(:,k);
end

%% Actually do the estimation
if 0
    dat=paramsPopulator('shadmuss');
    f=find(dat(:,4));
    lf=length(f);
    
    simAsm=repeatedSim(paramsPopulator(A(:,1:lf)),t,xvafsm,@armdynamics_inverted,0,1);
    simBsm=repeatedSim(paramsPopulator(B(:,1:lf)),t,xvafsm,@armdynamics_inverted,0,1);
    
    for k=1:length(f)
        k
        if k==1
            simABsm=repeatedSim(paramsPopulator(AB(:,1:lf,k)),t,xvafsm,@armdynamics_inverted,0,1);
        else
            simABsm(k,:)=repeatedSim(paramsPopulator(AB(:,1:lf,k)),t,xvafsm,@armdynamics_inverted,0,1);
        end
    end
    
    save('KICKsm_white.mat','simAsm','simBsm','simABsm')
end

if 1
    dat=paramsPopulator('burdet');
    f=find(dat(:,4));
    lf=length(f);
    
    simA=repeatedSim(paramsPopulator(A(:,1:lf)),t,xvaf,'reflex',0,1);
    simB=repeatedSim(paramsPopulator(B(:,1:lf)),t,xvaf,'reflex',0,1);
    
    for k=1:length(f)
        k
        if k==1
            simAB=repeatedSim(paramsPopulator(AB(:,1:lf,k)),t,xvaf,'reflex',0,1);
        else
            simAB(k,:)=repeatedSim(paramsPopulator(AB(:,1:lf,k)),t,xvaf,'reflex',0,1);
        end
    end
    save('KICK_white.mat','simA','simB','simAB')
end