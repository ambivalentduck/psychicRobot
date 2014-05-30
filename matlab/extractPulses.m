function extractPulses(k)

global tFit kpFit

figure(k)
clf
hold on

load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])

W(end,1)
lowerlim=0;
Wd(Wd<lowerlim)=lowerlim;

set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

dcats=[trials.disturbcat];
F=find((dcats>0)&(dcats<5)); %Just like title implies, only the pulses

avemeK=zeros(length(F),200);
avemeS=zeros(length(F),200);

for c=1:length(F)
    c/length(F)

    kk=F(c);
    onset=find(vecmag(trials(kk).v)>.05,1,'first');
    start=max(onset-35,1);
    onset2=find(vecmag(trials(kk+1).v)>.1,1,'first');
    xvaf1=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
    xvaf2=[trials(kk+1).x trials(kk+1).v trials(kk+1).a trials(kk+1).f];
    
    t=[trials(kk).t(start:end); trials(kk+1).t(1:onset2)]';
    t=t(1:min(length(t),trialInfo(kk).forceinds(1)+201));
    xvaf=[xvaf1(start:end,:); xvaf2(1:onset2,:)];
    
    tFit=trials(kk).t(trialInfo(kk).forceinds(1))+(0:.005:.005*length(Wd));
    kpFit=Wd;
    yK=extract(t,xvaf,'reflex');
    trials(kk).yK=yK;
    trials(kk).ty=t;
    
    tFit=1;
    kpFit=.8*W(end,1);
    yS=extract(t,xvaf,'reflex');
    trials(kk).yS=yS;
    
%     first=trialInfo(kk).forceinds(1);
%     tp=[trials(kk).t(first:end); trials(kk+1).t(1:onset2)]';
%     xvafp=[xvaf1(first:end,:); xvaf2(1:onset2,:)];
%     yp=extract(tp,xvafp,'reflex');
%     trials(kk).yp=yp;
%     trials(kk).typ=tp;

%     yoff=-xvaf(1,2)+trialInfo(kk).startcat;
%     xoff=-xvaf(1,1)+trialInfo(kk).endcat;
%     plot(xvaf1(:,1)+xoff,xvaf1(:,2)+yoff,'b')
%     plot(y(:,1)+xoff,y(:,2)+yoff,'r')
    si=sign(xvaf1(trialInfo(kk).forceinds(1)+40,2)-.5);
    toff=trials(kk).t(trialInfo(kk).forceinds(1));
    plot(t-toff,si*(xvaf(1:size(yS,1),2)-.5),'b')
    plot(t-toff,si*(yK(:,2)-.5),'r')
    plot(t-toff,si*(yS(:,2)-.5),'color',[245 208 76]/255)
    %plot(yp(:,1)+xoff,yp(:,2)+yoff,'g')
    plot([0 1],.008+[0 0],'k','linewidth',3)
    plot([0 1],-.012+[0 0],'k','linewidth',3)
    blah=find(((t-toff)>=0)&((t-toff)<=1));
    avemeS(c,1:length(blah))=si*(yS(blah,2)'-.5);
    avemeK(c,1:length(blah))=si*(yK(blah,2)'-.5);
    try
        delete(aveh)
    end
    aveh=[plot(linspace(0,1,length(blah)),mean(avemeS(1:c,1:length(blah))),'g','linewidth',3);
        plot(linspace(0,1,length(blah)),mean(avemeK(1:c,1:length(blah))),'c','linewidth',3)];
    xlim([0 1])
    drawnow
end
save(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'],'trials')

cleanup