clc
clear all

global kp0gain kp1gain

for k=1:4
    figure(k)
    clf
    hold on
    
    load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    kp0gain=W(end,1);
    kp1gain=W(end,1);
    set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)
    
    F=find(dcats>0);
    for c=1:length(F)
        c/length(F)
        
        kk=F(c);
        onset=find(vecmag(trials(kk).v)>.1,1,'first');
        start=max(onset-10,1);
        onset2=find(vecmag(trials(kk+1).v)>.1,1,'first');
        xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
        xvaf2=[trials(kk+1).x trials(kk+1).v trials(kk+1).a trials(kk+1).f];
        t=[trials(kk).t(start:end); trials(kk+1).t(1:onset2)];
        t=t'-t(1);
        trials(kk).tcat=t;
        xvaf=[xvaf(start:end,:); xvaf2(1:onset2,:)];

        y=extract(t,xvaf,'reflex');
        trials(kk).y=y;
        
        yoff=-xvaf(1,2)+starts(kk);
        xoff=-xvaf(1,1)+dcats(kk);
        plot(xvaf(:,1)+xoff,xvaf(:,2)+yoff,'b')
        plot(y(:,1)+xoff,y(:,2)+yoff,'r')
        drawnow
    end
    save(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'],'trials')
end