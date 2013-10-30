clc
clear all
close all

global kpgain
kpgain=1;
gray=.7;
vecmag_scale=.125;

for k=5 %1:3 %8
    name=num2str(k);
    if exist(['../Data/curlkick',name,'.mat'],'file')
        load(['../Data/curlkick',name,'.mat']);
        params=getSubjectParamsCurl(name);
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

    targlocs=zeros(max(targetcat),2);
    for c=1:max(targetcat)
        f=find(targetcat==c);
        x=zeros(length(f),2);
        for cc=1:length(f)
            x(cc,:)=trials(f(cc)).x(end,:);
        end
        targlocs(c,:)=mean(x);
    end
    
    for c=1:max(targetcat)
        if c==center
            rotmat{c}=[1 0;0 1];
        else
            angle=atan2(targlocs(c,2)-targlocs(center,2),targlocs(c,1)-targlocs(center,1))
            co=cos(angle);
            si=sin(angle);
            rotmat{c}=[co si;-si co];
            urotmat{c}=[co -si;si co];
        end
    end
    

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
    fclean=find(clean&(targetcat~=center));
    r=randperm(length(fclean));

    for c=1:60 %length(f)
        kk=fclean(r(c));

        t=trials(kk).t';
        t=t-t(1);

        onset=find(vecmag(trials(kk).v)>.1,1,'first');
        start=max(onset-10,1);

        xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
        
        x=trials(kk).x;
        dist=sum(vecmag(x(2:end,:)-x(1:end-1,:)));
        if dist>.25
            continue
        end

        subplot(3,2,1)
        hold on
        plot(xvaf(:,1),xvaf(:,2),'Color',[gray gray gray])
        %quiver(xvaf(1:10:end,1),xvaf(1:10:end,2),xvaf(1:10:end,7),xvaf(1:10:end,8),'b')
        %plot(y(:,1),y(:,2),'Color',[1 gray gray])
        axis equal

        subplot(3,2,2)
        hold on
        plot(xvaf(:,1),xvaf(:,2),'Color',[gray gray gray])
        %quiver(xvaf(1:10:end,1),xvaf(1:10:end,2),xvaf(1:10:end,7),xvaf(1:10:end,8),'b')
        %plot(y(:,1),y(:,2),'Color',[1 gray gray])
        axis equal

        vmy=vecmag(xvaf(:,3:4));
        xoff=xvaf(1,1:2);
        x=[xvaf(:,1)-xoff(1) xvaf(:,2)-xoff(2)];
        urx=x*urotmat{targetcat(kk)};
        rc=[urx(:,1) vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rcn=[urx(:,1) -vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rc=[rc(:,1)+xoff(1) rc(:,2)+xoff(2)];
        rcn=[rcn(:,1)+xoff(1) rcn(:,2)+xoff(2)];
        tpos=targlocs(targetcat(kk),:);
        cpos=targlocs(center,:);
                
        subplot(3,2,3)
        hold on
        plot([tpos(1),cpos(1)],[tpos(2), cpos(2)],'w')
        plot(rc(:,1),rc(:,2),'Color',[gray gray gray])
        plot(rcn(:,1),rcn(:,2),'Color',[gray gray gray])
        axis equal
        
        subplot(3,2,4)
        hold on
        plot([tpos(1),cpos(1)],[tpos(2), cpos(2)],'w')
        plot(rc(:,1),rc(:,2),'Color',[gray gray gray])
        plot(rcn(:,1),rcn(:,2),'Color',[gray gray gray])
        axis equal
        
        subplot(3,2,5)
        hold on
        plot(t-t(start),vmy,'Color',[gray gray gray])

        subplot(3,2,6)
        hold on
        plot(t-t(start),vmy,'Color',[gray gray gray])
    end

    problemtrials=[];
    for c=1:length(f)
        c/length(f)
        kk=f(c);

        onset=find(vecmag(trials(kk).v)>.1,1,'first');
        start=max(onset-10,1);

        onset2=find(vecmag(trials(kk+1).v)>.25,1,'first');

        xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];
        xvaf2=[trials(kk+1).x trials(kk+1).v trials(kk+1).a trials(kk+1).f];

        t=[trials(kk).t(start:end); trials(kk+1).t(1:onset2)];
        t=t'-t(1);

        xvaf=[xvaf(start:end,:); xvaf2(1:onset2,:)];

        mark=length(trials(kk).t(start:end));
        y=extract(t,xvaf,'reflex');
        trials(kk).y=y;
        lumps=findLumps(t,y,mark-10,0);
        %lumps=findLumps(t,y,mark-10,length(lumps));
        trials(kk).lumps=lumps;
        trials(kk).nlumps=length(lumps);
        if ((length(lumps)/t(end)>6)&&(length(lumps)>9))||~isfield(lumps(end),'ownership')
            [c kk]
            problemtrials(end+1,:)=[c kk];
            trials(kk).isproblem=1;
            continue
        end
        trials(kk).isproblem=0;
        
        [C,owned]=lumps2rgbk(lumps); 
        owned=find(owned)';
        owned=owned(owned<mark);


        figure(k)
        subplot(3,2,disturbcat(kk))
        plot(xvaf(1:mark,1),xvaf(1:mark,2),'k')
        %quiver(xvaf(1:10:end,1),xvaf(1:10:end,2),xvaf(1:10:end,7),xvaf(1:10:end,8),'b')
        for cc=owned
            plot(y(cc:cc+1,1),y(cc:cc+1,2),'-','Color',C(cc,:))
        end

        subplot(3,2,2+disturbcat(kk))
        vmy=vecmag(y(:,3:4));
        urx=[y(:,1)-y(1,1) y(:,2)-y(1,2)]*urotmat{targetcat(kk)};
        rc=[urx(:,1) -vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rc=[rc(:,1)+y(1,1) rc(:,2)+y(1,2)];
        for cc=owned
            plot(rc(cc:cc+1,1),rc(cc:cc+1,2),'-','Color',C(cc,:))
        end

        vmy=vecmag(xvaf(:,3:4));
        urx=[xvaf(:,1)-xvaf(1,1) xvaf(:,2)-xvaf(1,2)]*urotmat{targetcat(kk)};
        rc=[urx(:,1) vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rc=[rc(:,1)+xvaf(1,1) rc(:,2)+xvaf(1,2)];
        plot(rc(:,1),rc(:,2),'k')

        subplot(3,2,4+disturbcat(kk))
        plot(t(1:mark),vecmag(xvaf(1:mark,3:4)),'k')
        %plot(t(locs)-t(start),vecmag(xvaf(locs,3:4)),'rx')
        for cc=owned
            plot(t(cc:cc+1)',vecmag(y(cc:cc+1,3:4)),'-','Color',C(cc,:))
        end
    end
    subplot(3,2,1)
    axis off
    plot([-.05 .05],[.35 .35],'k')
    text(-.015,.337,'10cm')
    subplot(3,2,2)
    axis off
    subplot(3,2,3)
    axis off
    subplot(3,2,4)
    axis off
    subplot(3,2,5)
    axis off
    plot([0 0],[.4 .5],'k')
    text(-.05,.38,'10 cm/sec','Rotation',90)
    plot([0 1],[0 0],'k')
    text(.4,-.03,'1 sec')
    subplot(3,2,6)
    axis off
    backgray=.6;
    set(gcf,'Color',backgray*[1 1 1])
    drawnow
    cleanup
    
    nlumps=[trials.nlumps];
    bins=1:max(nlumps)
    h=hist(nlumps,bins);
    figure(k+10)
    bar(bins,h/sum(h))
    title(['Histogram of subunit count. Sub #',num2str(k)])
    ylabel('Relative Frequency')
end



