clc
clear all
close all

global kpgain
kpgain=1;
gray=.7;

for k=1 %:8
    name=num2str(k);
    if exist(['../Data/curlkick',name,'.mat'],'file')
        load(['../Data/curlkick',name,'.mat'])
    else
        input=load(['../Data/Data/input',name,'.dat']);
        output=load(['../Data/Data/output',name,'.dat']);
        traw=output(:,2);
        xvafraw=output(:,3:10);
        trial=output(:,1);
        f=find(trial==810,1,'last');

        params=getSubjectParamsCurl(name)

        targs=(min(input(:,3))+input(:,3))*pi+(min(input(:,4))+input(:,4));
        [blah1,blah2,targetcat]=unique(targs);
        [blah1,blah2,disturbcat]=unique(input(:,2));
        disturbcat(disturbcat==2)=0;
        disturbcat(disturbcat==3)=2;

        inds=1:f;
        trials=generalAddSubject(name,traw(inds),xvafraw(inds,:),trial(inds),targetcat,disturbcat,params);
        save(['../Data/curlkick',name,'.mat'],'trials','params')
    end

    targetcat=[trials.targetcat];
    disturbcat=[trials.disturbcat];

    set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

    h=hist(targetcat,unique(targetcat));
    [trash,center]=max(h);


    k
    figure(k)
    clf

    f=find((targetcat~=center)&(disturbcat~=0));

    nclean=4;
    clean=zeros(size(disturbcat));
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) disturbcat(1:end-(kk-1))];
    end
    clean=clean==0;
    fclean=find(clean);
    r=randperm(length(fclean));

    for c=1:60 %length(f)
        kk=fclean(r(c));

        t=trials(kk).t';
        t=t-t(1);

        onset=find(vecmag(trials(kk).v)>.1,1,'first');
        start=max(onset-10,1);

        xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
        %xvaf(:,8)=-xvaf(:,8);
        %y=extract(t(start:end),xvaf(start:end,:),'reflex');
        %trials(kk).y=y;

        subplot(2,2,1)
        hold on
        plot(xvaf(:,1),xvaf(:,2),'Color',[gray gray gray])
        %quiver(xvaf(1:10:end,1),xvaf(1:10:end,2),xvaf(1:10:end,7),xvaf(1:10:end,8),'b')
        %plot(y(:,1),y(:,2),'Color',[1 gray gray])
        axis equal
        
        subplot(2,2,2)
        hold on
        plot(xvaf(:,1),xvaf(:,2),'Color',[gray gray gray])
        %quiver(xvaf(1:10:end,1),xvaf(1:10:end,2),xvaf(1:10:end,7),xvaf(1:10:end,8),'b')
        %plot(y(:,1),y(:,2),'Color',[1 gray gray])
        axis equal

        subplot(2,2,3)
        hold on
        plot(t-t(start),vecmag(xvaf(:,3:4)),'Color',[gray gray gray])
        %plot(t(start:end)-t(start),vecmag(y(:,3:4)),'Color',[1 gray gray])
        
        subplot(2,2,4)
        hold on
        plot(t-t(start),vecmag(xvaf(:,3:4)),'Color',[gray gray gray])
        %plot(t(start:end)-t(start),vecmag(y(:,3:4)),'Color',[1 gray gray])
    end

    for c=1:length(f)
        kk=f(c);
        subplot(2,2,disturbcat(kk))
        hold on

        t=trials(kk).t';
        t=t-t(1);

        onset=find(vecmag(trials(kk).v)>.1,1,'first');
        start=max(onset-10,1);

        xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
        %xvaf(:,8)=-xvaf(:,8);
        y=extract(t(start:end),xvaf(start:end,:),'reflex');
        trials(kk).y=y;
        plot(xvaf(:,1),xvaf(:,2),'k')
        %quiver(xvaf(1:10:end,1),xvaf(1:10:end,2),xvaf(1:10:end,7),xvaf(1:10:end,8),'b')
        plot(y(:,1),y(:,2),'r')
        axis equal

        subplot(2,2,2+disturbcat(kk))
        hold on
        plot(t-t(start),vecmag(xvaf(:,3:4)),'k')
        %plot(t(locs)-t(start),vecmag(xvaf(locs,3:4)),'rx')
        plot(t(start:end)-t(start),vecmag(y(:,3:4)),'r')
    end
    subplot(2,2,1)
    axis off
    plot([-.05 .05],[.35 .35],'k')
    text(-.015,.337,'10cm')
    subplot(2,2,2)
    axis off
    subplot(2,2,3)
    axis off
    plot([0 0],[.4 .5],'k')
    text(-.05,.38,'10 cm/sec','Rotation',90)
    plot([0 1],[0 0],'k')
    text(.4,-.03,'1 sec')
    subplot(2,2,4)
    axis off

    drawnow
    cleanup
end


