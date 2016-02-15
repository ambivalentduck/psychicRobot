function extractReaches(k)

global kpgain

figure(k)
clf
hold on
axis equal

load(['../Data/curlkick/curlkick',num2str(k),'.mat'])
c=1;
mass=50;
kpgain=c;
[l1,l2,x0]=getSubjectParams(k);


set2dGlobals(l1,l2,x0,mass)

dcats=[trials.disturbcat];
tcats=[trials.targetcat];
unique(dcats)
F=find((dcats>0)&(tcats~=0));


for c=1:length(F)
    c/length(F)

    kk=F(c);
    start=curlKickOnsetDetector(trials(kk));
    onset0=curlKickOnsetDetector(trials(kk-1));
    onset2=curlKickOnsetDetector(trials(kk+1))+30;

    % The extraction is sensitive to the assumption of x0, but it quickly settles
    xvaf0=[trials(kk-1).x trials(kk-1).v trials(kk-1).a trials(kk-1).f];
    xvaf1=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
    xvaf2=[trials(kk+1).x trials(kk+1).v trials(kk+1).a trials(kk+1).f];
    
    t=[trials(kk-1).t(onset0:end); trials(kk).t; trials(kk+1).t(1:onset2)]';
    t=t-t(start);
    xvaf=[xvaf0(onset0:end,:); xvaf1; xvaf2(1:onset2,:)];
    
    y=extract(t,xvaf,'reflex');
    trials(kk).y=y;
    trials(kk).ty=t;
    trials(kk).i0=length(trials(kk-1).t)-onset0+1;
    trials(kk).if=length(trials(kk-1).t)-onset0+length(trials(kk).t);
    
    if 1
    yoff=trials(kk).disturbcat/2;
    xoff=0;%-xvaf(1,1); %+trials(kk).targetcat/3;
    plot(xvaf1(:,1)+xoff,xvaf1(:,2)+yoff,'b')
    inds=trials(kk).i0:trials(kk).if;
    plot(y(inds,1)+xoff,y(inds,2)+yoff,'r')
    drawnow
    end
end
save(['../Data/curlkick/curlkick',num2str(k),'Y.mat'],'trials')

cleanup