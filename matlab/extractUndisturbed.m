function extractUndisturbed(SUBNUM)

global kp0gain kp1gain

load(['../Data/Data_pulse/pulse',num2str(SUBNUM),'.mat'])
load(['../Data/Data_pulse/pulse',num2str(SUBNUM),'W.mat'])

kp0gain=W(end,1);
kp1gain=W(end,1);
set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

starts=[trialInfo.startcat];
ends=[trialInfo.endcat];
clean=[trialInfo.clean];

figure(SUBNUM+400)
clf
hold on

f=find(clean);

for k=1:length(f)
    k/length(f)
    kk=f(k);
    if kk==1
        continue
    end
        
    onset=find(vecmag(trials(kk).v)>.05,1,'first');
    start=max(onset-35,1);
    xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
    undisturbed(kk).yinds=start:length(trials(kk).t);
    t=trials(kk).t(start:end);
    undisturbed(kk).ty=t;
    t=t'-t(1);
    y=extract(t,xvaf,'reflex');
    undisturbed(kk).y=y;
    x=trials(kk).x;

    [bx,trash1,trash2,trash2,statsx]=regress(x(:,2),[ones(size(x,1),1) x(:,1)]);
    undisturbed(kk).R2x=statsx(1);
    mx=[min(x(:,1)) max(x(:,1))];

    [by,trash1,trash2,trash2,statsy]=regress(y(:,2),[ones(size(y,1),1) y(:,1)]);
    undisturbed(kk).R2y=statsy(1);
    my=[min(y(:,1)) max(y(:,1))];

    yoff=trialInfo(kk).reachcat*.05;
    plot(x(:,1),x(:,2)+yoff,'b.')
    plot(y(:,1),y(:,2)+yoff,'r.')
    plot(mx,bx(2)*mx+bx(1)+yoff,'b-')
    plot(my,by(2)*my+by(1)+yoff,'r-')
end
axis equal

save(['../Data/Data_pulse/pulse',num2str(SUBNUM),'U.mat'],'undisturbed')