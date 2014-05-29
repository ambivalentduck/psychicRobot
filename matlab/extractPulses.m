function extractPulses(k)

global kp0gain kp1gain

figure(k)
clf
hold on

scale=.8;

load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
kp0gain=scale*W(end,1);
kp1gain=scale*W(end,1);
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
    
    first=trialInfo(kk).forceinds(1);
    tp=[trials(kk).t(first:end); trials(kk+1).t(1:onset2)]';
    xvafp=[xvaf1(first:end,:); xvaf2(1:onset2,:)];
    yp=extract(tp,xvafp,'reflex');
    trials(kk).yp=yp;
    trials(kk).typ=tp;

    yoff=-xvaf(1,2)+trialInfo(kk).startcat;
    xoff=-xvaf(1,1)+trialInfo(kk).endcat;
    plot(xvaf1(:,1)+xoff,xvaf1(:,2)+yoff,'b')
    plot(y(:,1)+xoff,y(:,2)+yoff,'r')
    plot(yp(:,1)+xoff,yp(:,2)+yoff,'g')
    drawnow
end
save(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'],'trials')

cleanup