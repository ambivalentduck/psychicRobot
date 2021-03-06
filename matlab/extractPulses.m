function extractPulses(k)

global kpgain massgain reflexcontrib

massgain=1;
 
kg=[.18];
rc=[.07 .1 .1 .1 .1 .1 .1 .1];
kpgain=kg(k);
reflexcontrib=rc(k);

figure(k)
clf
hold on
axis equal

load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])

set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

dcats=[trials.disturbcat];
F=find((dcats>0)&(dcats<5)); %Just like title implies, only the pulses

for c=1:length(F)
    c/length(F)

    kk=F(c);
    start=onsetDetector(trials(kk));
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