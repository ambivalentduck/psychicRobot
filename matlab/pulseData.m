clc
clear all
close all

global kpgain
kpgain=1;
gray=.7;
vecmag_scale=.125;
WIDTH=3;

subnums=[324 789 300 301];
kpgains=[.4 .25 1 .5];

for k=[1 2 3 4] %:length(subnums)
    if ~exist(['../Data/Data_pulse/pulse',num2str(k),'.mat'],'file')||1
        name=num2str(subnums(k));
        input=load(['../Data/Data_pulse/input',name,'.dat']);
        output=load(['../Data/Data_pulse/output',name,'.dat']);
        params=getSubjectParams(name);
        traw=output(:,2);
        xvafraw=output(:,3:10);
        trial=output(:,1);
        lastTrial=max(trial)-1;
        f=find(trial==lastTrial,1,'last');
        samprate=1/mean(gradient(traw))
        if samprate<300 %our filtering method apparently has nonlinear algorithmic complexity and goes from seconds to hours+.
            inds=1:f;
        else
            inds=1:5:f;
        end
        trials=generalAddSubject(name,traw(inds),xvafraw(inds,:),trial(inds),params);
        figure(10)
        clf
        subplot(2,1,1)
        hold on
        TN=18;
        plot(trials(TN).x(:,1),trials(TN).x(:,2))
        quiver(trials(TN).x(:,1),trials(TN).x(:,2),trials(TN).f(:,1),trials(TN).f(:,2))
        axis equal
        subplot(2,1,2)
        plot(trials(TN).t,vecmag(trials(TN).f))

        %Because trial 1 starts from an unknown location/where the robot is
        %turned on
        trials(1).dist=.15;
        trials(1).dir=-1;
        trials(1).disturbcat=0;

        cats=[1 0 2];
        dists=[.15 .3];

        for c=2:length(trials)
            [trash,ind]=min(abs(norm(trials(c).x(end,:)-trials(c).x(1,:))-dists));
            trials(c).dist=dists(ind);
            trials(c).dir=sign(trials(c).x(end,1)-trials(c).x(1,1));
            trials(c).disturbance=input(c,[4 5 6]);
            v=find(input(c,[4 5 6]));
            if isempty(v)
                trials(c).disturbcat=0;
            else
                trials(c).disturbcat=2*v-(input(c,3+v)>=0);
            end
        end
        unique([trials.dist])
        unique([trials.dir])
        unique([trials.disturbcat])
        save(['../Data/Data_pulse/pulse',num2str(k),'.mat'],'trials','params')
    end
end

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'.mat']);

    udists=unique([trials.dist]);
    udirs=unique([trials.dir]);
    sumudists=sum(udists);
    disturbcat=[trials.disturbcat];
    udisturbcat=unique(disturbcat)

    set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

    k
    figure(k)
    clf
    UD=-.3;
    LR=1.8*sumudists;
    Toff=[0 1.9];
    
    hold on
    for dist=udists
        for dir=udirs
            distn=find(dist==udists);
            dirn=find(dir==udirs);
            xoff=LR(distn,dirn)
            yoff=0;
            for cc=0:2
                plot(xoff,yoff+cc*UD,'wo','markerfacecolor','w')
                plot([0,dir*dist]+xoff,[0 0]+yoff+cc*UD,'w','LineWidth',3)
            end
        end
    end
    axis equal
    axis off

    nclean=4;
    clean=zeros(size(disturbcat));
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) disturbcat(1:end-(kk-1))];
    end
    clean=clean==0;
    fclean=find(clean);
    r=randperm(length(fclean));

    for c=1:90
        kk=fclean(r(c));
        if kk==1
            continue
        end

        t=trials(kk).t';
        t=t-t(1);

        onset=find(vecmag(trials(kk).v)>.1,1,'first');
        start=max(onset-10,1);

        xvaf=[trials(kk).x trials(kk).v trials(kk).a trials(kk).f];

        x=trials(kk).x;
        dist=sum(vecmag(x(2:end,:)-x(1:end-1,:)));
        if dist>1.2*trials(kk).dist %Throw away baseline examples where they sneezed or something.
            continue
        end

        distn=find(trials(kk).dist==udists);
        dirn=find(trials(kk).dir==udirs);
        xoff=xoffsets(distn,dirn);
        yoff=0;

        for cc=yoffsets
            plot(xvaf(:,1)-xvaf(1,1)+xoff,xvaf(:,2)-xvaf(1,2)+yoff+cc,'Color',[gray gray gray])
        end
    end
    drawnow
    colscheme=colorScheme(10);
    problemtrials=[];
    anecdotes=zeros(4,2);

    f=find([trials.disturbcat]);

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
        %[y,lumps,breakp(c),bestKp(c,:)]=stiffnessExtract(t,xvaf,0);
        lumps=findLumps(t,y,mark-10,0);
        trials(kk).y=y;
        trials(kk).cumdist=[0; cumsum(vecmag(trials(kk).y(2:end,1:2)-trials(kk).y(1:end-1,1:2)))];
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

        distn=find(trials(kk).dist==udists);
        dirn=find(trials(kk).dir==udirs);
        xoff=xoffsets(distn,dirn);
        yoff=yoffsets(trials(kk).disturbcat);
        TOF=Toff(distn);

        x=[xvaf(:,1)-xvaf(1,1)+xoff xvaf(:,2)-xvaf(1,2)+yoff]; %New plotting needs to start at 0.

        %X
        plot(x(1:mark,1),x(1:mark,2),'k')

        %f
        qscale=200;
        qskip=5;
        quiver(x(1:qskip:end,1),x(1:qskip:end,2),xvaf(1:qskip:end,7)/qscale,xvaf(1:qskip:end,8)/qscale,0,'Color',[.54 .27 .075])

        %Y
        yplot=[y(:,1)-y(1,1)+xoff y(:,2)-y(1,2)+yoff];
        for cc=owned
            plot(yplot(cc:cc+1,1),yplot(cc:cc+1,2),'-','Color',C(cc,:))
        end

        drawnow
    end
    save(['../Data/pulse',num2str(k),'g.mat'],'trials');

    set(gcf,'Position',[66 1 1215 945]);

    %30 cm
    plot([-.3 0]+xoffsets(2,1),.02*[1 1],'w','LineWidth',3)
    text(-.17+xoffsets(2,1),.03,'30 cm','color',[1 1 1])
    %10 N
    plot(.01*[1 1],UD+.05+(1/200)*[0 10],'w','LineWidth',3)
    text(.015,UD+.075,'10 N','color',[1 1 1])

    backgray=.6;
    set(gcf,'Color',backgray*[1 1 1])
    drawnow

    gf=getframe(gcf);
    imwrite(gf.cdata,['./posterprep/summary',num2str(k),'pulse.png'])
end

%At the end of the red, which trace "wins" on lowest deflection from
%average? Two questions...what's the average gray at the end of the
%red. And what are the deflections at the end of the red.

%for each disturbed, find a gray of the same dir and dist, and grab its
%deflect there.

deflects={};
deflectshand={};
basedeflects={};
groups={};
signcat=[-1 1];
for k=1:4
    load(['../Data/pulse',num2str(k),'g.mat']);
    disturbcat=[trials.disturbcat];
    f=find(disturbcat);
    dirs=[trials.dir];
    dists=[trials.dist];
    nlumps=[trials.nlumps];
    nclean=4;
    clean=zeros(size(dirs));
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) disturbcat(1:end-(kk-1))];
    end
    clean=clean==0;
    for L=1:WIDTH
        if L>nlumps(c)
            continue
        end

        groups{k,L}=k*ones(size(f));
        deflects{k,L}=f;
        deflectshand{k,L}=f;
        basedeflects{k,L}=f;
        for c=1:length(f)
            kk=f(c);
            onset=find(vecmag(trials(kk).v)>.1,1,'first');
            start=max(onset-10,1);
            ind=trials(kk).lumps(L).inds(end);
            deflects{k,L}(c)=signcat(disturbcat(kk))*(trials(kk).y(ind,2)-.5);
            deflectshand{k,L}(c)=signcat(disturbcat(kk))*(trials(kk).x(ind+start,2)-.5);
            f2=find(clean&(dists==trials(kk).dist)&(dirs==trials(kk).dir));
            shuf=randperm(length(f2));
            [trash,ind2]=min(abs(trials(f2(shuf(1))).x(:,1)-trials(k).x(ind,1)));
            basedeflects{k,L}(c)=trials(f2(shuf(1))).x(ind2,2)-.5;
        end
    end
end

pos=0:4;
pos=[pos-.2 pos pos+.2];

cols=zeros(15,3);
for k=1:5
    cols(k,:)=[0 0 0];
end
for k=6:10
    cols(k,:)=[gray gray gray];
end

colscheme=colorScheme(10);

for k=1:WIDTH
    for cc=11:15
        cols(cc,:)=colscheme(k,:);
    end
    subplot(NSP,WIDTH,WIDTH*(NSP-1)+k)
    groups2=[groups{:,k}]';
    groups2=[groups2; 0*groups2];
    groups2=[groups2; groups2+5;groups2+10];
    dh=[deflectshand{:,k}]';
    bd=[basedeflects{:,k}]';

    d=[deflects{:,k}]';
    dat2=[dh; dh; bd; bd; d; d];
    %boxplot(abs(dat2),groups2,'orientation','vertical','plotstyle','compact','positions',pos,'colors',cols)
    boxplot(abs(dat2),groups2,'orientation','vertical','notch','on','positions',pos,'colors',cols,'width',.15)
    set(gca,'Color',backgray*[1 1 1])
    set(gca,'xtick',0:4)
    set(gca,'xticklabel',{'All','1','2','3','4'})
    if k==1
        ylabel('|Deflection|')
    else
        ylabel('')
    end
    if k==2
        xlabel('Subject')
    end
end



subplot(htop)
apos=get(gca,'Position');
%set(gca,'Position',[0 0 1215 apos(4)*1215/apos(3)])
set(gca,'Position',[.05 1-.7 .9 .8])

gf=getframe(gcf);
imwrite(gf.cdata,['./posterprep/summary',num2str(k),'pulse.png'])