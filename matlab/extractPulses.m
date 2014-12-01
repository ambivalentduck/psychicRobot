function extractPulses(k)

global kpgain

figure(k)
clf
hold on
axis equal

scale=.8;

load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
params
%Subject 1
% params.mass=50;
% params.shoulder(2)=.48;
% params.l1=.28;
% params.l2=.3;

%Subject 2
kpgain=.16

%Subject 3 and 4
%kpgain=scale*W(end,1)

set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

offsetForce=onset;

dcats=[trials.disturbcat];
F=find((dcats>0)&(dcats<5)); %Just like title implies, only the pulses

for c=1:length(F)
    c/length(F)

    kk=F(c);
    onset=find(vecmag(trials(kk).v)>.05,1,'first');
    start=max(onset-35,1);
    onset2=find(vecmag(trials(kk+1).v)>.1,1,'first');
    xvaf1=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
    xvaf2=[trials(kk+1).x trials(kk+1).v trials(kk+1).a trials(kk+1).f];
    
    t=[trials(kk).t(start:end); trials(kk+1).t(1:onset2)]';
    xvaf=[xvaf1(start:end,:); xvaf2(1:onset2,:)];
    y=extract(t,xvaf,'reflex');
    trials(kk).y=y;
    trials(kk).ty=t;
    
    yoff=-xvaf(1,2)+trialInfo(kk).startcat/3;
    xoff=-xvaf(1,1)+trialInfo(kk).endcat/3;
    plot(xvaf1(:,1)+xoff,xvaf1(:,2)+yoff,'b')
    plot(y(:,1)+xoff,y(:,2)+yoff,'r')
    drawnow
end
save(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'],'trials')

cleanup