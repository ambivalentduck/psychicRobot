clc
clear all
close all

global kpgain
kpgain=1;
gray=.7;
vecmag_scale=.125;
starts={};
ends={};

for k=1:8
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
    NSP=3;
    UD=-.4;
    LR=.55;
    Toff=[0 1.9];

    subplot(NSP,1,1:NSP-1)
    hold on
    for cc=LR*(0:3)
        for ccc=1:4
            tpos=targlocs(ccc,:);
            cpos=targlocs(center,:);
            plot([tpos(1),cpos(1)]+cc,[tpos(2), cpos(2)]+UD,'w','LineWidth',3)
        end
    end
    axis off
    subplot(NSP,1,NSP)
    axis off
    hold on

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
        if dist>.25 %Throw away baseline examples where they sneezed or something.
            continue
        end

        subplot(NSP,1,1:NSP-1)
        hold on
        for cc=LR*(0:3)
            plot(xvaf(:,1)+cc,xvaf(:,2),'Color',[gray gray gray])
        end

        vmy=vecmag(xvaf(:,3:4));
        xoff=xvaf(1,1:2);
        x=[xvaf(:,1)-xoff(1) xvaf(:,2)-xoff(2)];
        urx=x*urotmat{targetcat(kk)};
        rc=[urx(:,1) vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rcn=[urx(:,1) -vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rc=[rc(:,1)+xoff(1) rc(:,2)+xoff(2)];
        rcn=[rcn(:,1)+xoff(1) rcn(:,2)+xoff(2)];


        for cc=LR*[1 3]
            plot(rc(:,1)+cc,rc(:,2)+UD,'Color',[gray gray gray])
            plot(rcn(:,1)+cc,rcn(:,2)+UD,'Color',[gray gray gray])
        end
        axis equal

        subplot(NSP,1,NSP)
        hold on
        for cc=1:2
            plot(t(start:end)-t(start)+Toff(cc),vmy(start:end),'Color',[gray gray gray])
        end
    end
    drawnow

    colscheme=colorScheme(10);
    problemtrials=[];
    anecdotes=zeros(4,2);

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
        trials(kk).tcat=t;
        xvaf=[xvaf(start:end,:); xvaf2(1:onset2,:)];

        mark=length(trials(kk).t(start:end));
        y=extract(t,xvaf,'reflex');
        trials(kk).y=y;
        trials(kk).cumdist=[0; cumsum(vecmag(trials(kk).y(2:end,1:2)-trials(kk).y(1:end-1,1:2)))];
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

        subplot(NSP,1,1:NSP-1)
        XOF=2*LR*(disturbcat(kk)-1);
        TOF=Toff(disturbcat(kk));

        %X, f, gray
        plot(xvaf(1:mark,1)+XOF,xvaf(1:mark,2),'k')
        qscale=200;
        qskip=5;
        %quiver(xvaf(1:qskip:end,1)+XOF,xvaf(1:qskip:end,2),xvaf(1:qskip:end,7)/qscale,xvaf(1:qskip:end,8)/qscale,0,'Color',[1 .46 .094])
        quiver(xvaf(1:qskip:end,1)+XOF,xvaf(1:qskip:end,2),xvaf(1:qskip:end,7)/qscale,xvaf(1:qskip:end,8)/qscale,0,'Color',[.54 .27 .075])

        %X, Y, gray
        plot(xvaf(1:mark,1)+LR+XOF,xvaf(1:mark,2),'k')
        for cc=owned
            plot(y(cc:cc+1,1)+LR+XOF,y(cc:cc+1,2),'-','Color',C(cc,:))
        end

        %Anecdotes
        if (sum(length(lumps)==[3 4])>0)&&~anecdotes(targetcat(kk),disturbcat(kk))
            yls=(y(1,1:2)'*ones(1,length(t)))';
            for cc=1:length(lumps)
                inds=lumps(cc).inds;
                yl=[lumps(cc).y(:,1)-lumps(cc).y(1,1) lumps(cc).y(:,2)-lumps(cc).y(1,2)];
                offset=y(inds(floor(end/2+.5)),1:2)-yl(floor(end/2+.5),1:2);
                yls(inds,:)=yls(inds,:)+yl;
                yls(inds(end)+1:end,:)=[yls(inds(end)+1:end,1)+yl(end,1) yls(inds(end)+1:end,2)+yl(end,2)];
                ylp=[yl(:,1)+offset(1) yl(:,2)+offset(2)];
                plot(ylp(:,1)+XOF,ylp(:,2)+UD,'-','Color',colscheme(cc,:))
                %plot(yls(inds,1),yls(inds,2),'.','Color',colscheme(cc,:))
            end
            for cc=owned
                plot(y(cc:cc+1,1)+XOF,y(cc:cc+1,2)+UD,'-','Color',C(cc,:),'Linewidth',2)
            end
            anecdotes(targetcat(kk),disturbcat(kk))=1;
        end

        %Clover-leaf
        vmy=vecmag(y(:,3:4));
        urx=[y(:,1)-y(1,1) y(:,2)-y(1,2)]*urotmat{targetcat(kk)};
        trials(kk).yurot=urx;
        rc=[urx(:,1) -vmy*vecmag_scale]*rotmat{targetcat(kk)};
        rc=[rc(:,1)+y(1,1) rc(:,2)+y(1,2)];
        for cc=owned
            plot(rc(cc:cc+1,1)+LR+XOF,rc(cc:cc+1,2)+UD,'-','Color',C(cc,:))
        end
        for cc=1:length(lumps)
            inds=lumps(cc).inds;
            yl=[lumps(cc).y(:,1)-lumps(cc).y(1,1) lumps(cc).y(:,2)-lumps(cc).y(1,2)];
            offset=y(inds(floor(end/2+.5)),1:2)-yl(floor(end/2+.5),1:2);
            ylp=[yl(:,1)+offset(1) yl(:,2)+offset(2)];
            ury=[ylp(:,1)-y(1,1) ylp(:,2)-y(1,2)]*urotmat{targetcat(kk)};
            vm=vecmag(lumps(cc).y(:,3:4));
            rc=[ury(:,1) vm*vecmag_scale]*rotmat{targetcat(kk)};
            rc=[rc(:,1)+y(1,1) rc(:,2)+y(1,2)];
            plot(rc(:,1)+LR+XOF,rc(:,2)+UD,'color',colscheme(cc,:))
        end
        %         vmy=vecmag(xvaf(:,3:4));
        %         urx=[xvaf(:,1)-xvaf(1,1) xvaf(:,2)-xvaf(1,2)]*urotmat{targetcat(kk)};
        %         rc=[urx(:,1) vmy*vecmag_scale]*rotmat{targetcat(kk)};
        %         rc=[rc(:,1)+xvaf(1,1) rc(:,2)+xvaf(1,2)];
        %         plot(rc(:,1)+LR+XOF,rc(:,2)+UD,'k')

        subplot(NSP,1,NSP)
        plot(t(1:mark)+TOF,vecmag(xvaf(1:mark,3:4)),'k')
        %plot(t(locs)-t(start),vecmag(xvaf(locs,3:4)),'rx')
        for cc=owned
            plot(t(cc:cc+1)'+TOF,vmy(cc:cc+1)-.6,'-','Color',C(cc,:))
        end
        drawnow
    end
    save(['../Data/curlkick',name,'g.mat'],'trials');

    set(gcf,'Position',[66 1 1215 945]);
    hbottom=subplot(NSP,1,NSP)
    htop=subplot(NSP,1,1:NSP-1)
    subplot(htop)
    apos=get(gca,'Position');
    %set(gca,'Position',[0 0 1215 apos(4)*1215/apos(3)])
    set(gca,'Position',[.05 1-.7 .9 .8])
    %10 cm
    plot([-.03 -.03],[.15 .25],'w','LineWidth',3)
    text(-.045,.165,'10 cm','color',[1 1 1],'Rotation',90)
    %10 N
    plot([0 .05],.67*[1 1],'w','LineWidth',3)
    text(0,.68,'10 N','color',[1 1 1])
    % 1/8 m/s
    plot([.55 .63],.28*[1 1],'w','LineWidth',3)
    text(.56,.293,'1 m/s','color',[1 1 1])
    
    subplot(hbottom)
    set(gca,'Position',[.1 .15 .8 .3])
    plot([0 1],[0 0]-.6,'w','LineWidth',3)
    plot([0 0],[0 .2]-.6,'w','LineWidth',3)
    text(.95,-.65,'1 s','color',[1 1 1])
    text(-.02,-.61,'20 cm/s','color',[1 1 1],'Rotation',90)
    plot([0 1]+Toff(2),[0 0]-.6,'w','LineWidth',3)
    plot([0 0]+Toff(2),[0 .2]-.6,'w','LineWidth',3)
    text(.95+Toff(2),-.65,'1 s','color',[1 1 1])
    text(-.02+Toff(2),-.61,'20 cm/s','color',[1 1 1],'Rotation',90)
    
    backgray=.6;
    set(gcf,'Color',backgray*[1 1 1])
    drawnow
    
    gf=getframe(gcf);
    imwrite(gf.cdata,['./posterprep/summary',num2str(k),'.png'])

    nlumps=[trials.nlumps];
    bins=1:max(nlumps)
    h=hist(nlumps,bins);
    figure(k+10)
    bar(bins,h/sum(h))
    title(['Histogram of subunit count. Sub #',num2str(k)])
    ylabel('Relative Frequency')
    
    targdists=vecmag([targlocs(:,1)-targlocs(center,1) targlocs(:,2)-targlocs(center,2)]);
    nl=[trials(f).nlumps];
    ip=[trials(f).isproblem];
    for c=1:5
        f2=find((nl>=c)&~ip);
        starts{k,c}=zeros(length(f2),3);
        ends{k,c}=zeros(length(f2),3);
        for cc=1:length(f2)
            kk=f(f2(cc));
            inds=trials(kk).lumps(c).inds;
            rat=targdists(targetcat(kk));
            starts{k,c}(cc,1:3)=[trials(kk).cumdist(inds(1))/rat trials(kk).yurot(inds(1),1)/rat trials(kk).tcat(inds(1))];
            ends{k,c}(cc,1:3)=[trials(kk).cumdist(inds(end))/rat trials(kk).yurot(inds(end),1)/rat trials(kk).tcat(inds(end))];
        end
    end
end

submax=8;



subnums={};
for k=1:8
    for c=1:5
        starts{k,c}=starts{k,c}';
        ends{k,c}=ends{k,c}';
        subnums{k,c}=k*ones(size(starts{k,c}));
    end
end

figure(666)
clf
subplot(5,3,14)
nlumps=zeros(length(nlumps),submax);
groups=nlumps;
for k=1:submax
    load(['../Data/curlkick',num2str(k),'g.mat']);
    nlumps(:,k)=[trials.nlumps]';
    groups(:,k)=ones(45,1)*k;
    names{k+1}=num2str(k);
end
names{1}='All';
groups2=[groups(:); 0*groups(:)];
nlumps2=[nlumps(:); nlumps(:)];
%boxplot(nlumps2,groups2,'orientation','horizontal','notch','on','labels',names)
boxplot(nlumps2,groups2,'orientation','vertical','notch','on','labels',names)
%boxplot(nlumps2,groups2,'orientation','vertical','plotstyle','compact','labels',names)
xlabel('Subject Number')
ylabel('Submovements Detected')

yl=[3 1.5 2.5];
tits={'Cumulative Distance Ratio','Progress Ratio','Time, s'};
pos=0:8;
pos=[pos-.1 pos+.1];

cols=zeros(18,3);
% for k=1:18
%     if mod(k,2)
%         cols(k,:)=[0 0 1];
%     else
%         cols(k,:)=[1 0 0];
%     end
% end
cols(1:9,3)=1;
cols(10:18,1)=1;

for k=1:4
    dats=[starts{:,k}]';
    date=[ends{:,k}]';
    sms=[subnums{:,k}]';
    for met=1:3
        subplot(5,3,met+3*(k-1))
        dat2=[dats(:,met);dats(:,met);date(:,met);date(:,met)];
        f=find((~imag(dat2))&(dat2<50));
        sms2=[sms(:,met);0*sms(:,met);sms(:,met)+9;0*sms(:,met)+9];
        boxplot(dat2(f),sms2(f),'orientation','vertical','plotstyle','compact','positions',pos,'colors',cols)
        if k==1
            title(tits{met})
        end
        if met==1
            ylabel(['Submovement ',num2str(k)])
        else
            ylabel('')
        end
        ylim([0 yl(met)])
        set(gca,'xtick',0:8)
        set(gca,'xticklabel',names)
    end
end
    
   




