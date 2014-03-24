clc
clear all

global kp0gain kp1gain

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])

    % Categorize by start/end pair
    figure(1)
    clf
    hold on
    clustermeS=zeros(length(trials)-1,1);
    clustermeE=clustermeS;
    for t=1:length(trials)-1
        plot(trials(t+1).x([1 end],1),trials(t+1).x([1 end],2),'.')
        clustermeS(t)=trials(t+1).x(1,1);
        clustermeE(t)=trials(t+1).x(end,1);
    end
    [cats,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
    means=sort(means);
    plot(means,.5*[1 1 1],'rx')
    axis equal

    starts=zeros(length(trials),1);
    ends=starts;
    for t=1:length(trials)
        [trash,starts(t)]=min(abs(trials(t).x(1,1)-means));
        [trash,ends(t)]=min(abs(trials(t).x(end,1)-means));
    end
    starts(1)=-inf; %never use this reach...

    reachcat=10*starts+ends;
    urc=unique(reachcat);
    urc=urc(2:end);

    dcats0=[trials.disturbcat];

    if 0
        warning off all
        %No clue why, but the input file for #3 was whack. Found it through
        %commit history
        for kk=1:length(trials)
            if sum(vecmag(trials(kk).f))<30
                trials(kk).disturbance=[0 0 0];
            else
                [val,ind]=findpeaks(abs(trials(kk).f(:,1)),'minpeakheight',3);
                if length(ind)>2
                    trials(kk).disturbance=[0 0 1.5];
                else
                    [val,ind]=findpeaks(abs(trials(kk).f(:,2)),'minpeakheight',8,'npeaks',1);
                    if length(ind)>=1
                        if (sign(trials(kk).f(ind,2))*sign(trials(kk).x(end,1)-trials(kk).x(1,1)))<0
                            mag=15;
                        else
                            mag=-15;
                        end
                        if ((trials(kk).x(ind-28,1)-trials(kk).x(1,1))/(trials(kk).x(end,1)-trials(kk).x(1,1)))>.3
                            trials(kk).disturbance=[0 mag 0];
                        else
                            trials(kk).disturbance=[mag 0 0];
                        end
                    else

                    end
                end
            end
            v=find(trials(kk).disturbance);
            if isempty(v)
                trials(kk).disturbcat=0;
            else
                trials(kk).disturbcat=2*v-(trials(kk).disturbance(v)>=0);
            end
        end
        warning on all

        figure(500)
        clf
        plot(dcats0,dcats+.3*rand(size(dcats)),'.')
        return
    end

    % Use unperturbed examples for a reference trajectory for each start/end pair.
    dcats=[trials.disturbcat];

    nclean=3;
    clean=0*dcats;
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
    end
    clean=(clean==0);

    figure(2)
    clf
    hold on
    yoffset=.1;
    for U=1:length(urc)
        f=find(clean'&(reachcat==(urc(U))));
        rp=randperm(length(f));
        for krp=1:length(rp)
            fk=rp(krp);
            if krp<=20
                catme(U,krp).x=[trials(f(fk)).x(:,1) trials(f(fk)).v(:,1) trials(f(fk)).a(:,1)]';
                catme(U,krp).y=[trials(f(fk)).x(:,2) trials(f(fk)).v(:,2) trials(f(fk)).a(:,2)]';
            end
            plot(trials(f(fk)).x(:,1),trials(f(fk)).x(:,2)+yoffset*U,'linewidth',.01)
            plot(trials(f(fk)).x(end,1),trials(f(fk)).x(end,2)+yoffset*U,'rx')
        end
        st=floor(urc(U)/10);
        en=mod(urc(U),10);
        catx=[catme(U,:).x];
        caty=[catme(U,:).y];
        X=linspace(means(st),means(en));
        Y=gaussianWeightedRegression(catx(1,:)',caty',X,1000);
        plot(X,Y(:,1)+yoffset*U,'g')
    end
    axis equal

    yoffset=.3;
    figure(3)
    clf
    hold on
    for U=1:length(urc)
        st=floor(urc(U)/10);
        en=mod(urc(U),10);
        catx=[catme(U,:).x];
        caty=[catme(U,:).y];
        X=linspace(means(st),means(en));
        Y=gaussianWeightedRegression(catx(1,:)',caty',X,1000);
        plot(catx(1,:),caty(2,:)+yoffset*U,'b.','markersize',.01)
        plot(X,Y(:,2)+yoffset*U,'g')
    end

    % Find onset of forces
    figure(4)
    clf
    hold on
    k__=0;
    onsets=zeros(sum((dcats>0)&(dcats<5)),1);
    for U=1:length(urc)
        f=find((dcats>0)&(dcats<5)&(reachcat'==(urc(U))));
        for fk=1:length(f)
            [val,ind]=findpeaks(abs(trials(f(fk)).f(:,2)),'minpeakheight',5,'npeaks',1);
            if isempty(ind)
                continue
                disp(['Lowering the min for ',num2str(f(fk))])
                [val,ind]=findpeaks(abs(trials(f(fk)).f(:,2)),'minpeakheight',0,'npeaks',1);
            end
            time=trials(f(fk)).t-trials(f(fk)).t(ind);
            rforce=sign(trials(f(fk)).f(ind,2))*trials(f(fk)).f(:,2);
            ft=find(time<-.2);
            blforce=mean(rforce(ft));
            plot(time,blforce*ones(size(time))+20*U,'r');
            onset=find(rforce(1:ind)<=blforce,1,'last');
            plot(time,rforce+20*U,'linewidth',.01)
            plot(time(onset),rforce(onset)+20*U,'g.')
            k__=k__+1;
            onsets(k__)=time(onset);
        end
    end
    onset=mean(onsets)

    set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

    figure(5)
    clf
    hold on
    fitlower=0;
    fitupper=.25;
    for U=1:length(urc)
        f=find((dcats>0)&(dcats<5)&(reachcat'==(urc(U))));
        catx=[catme(U,:).x]';
        caty=[catme(U,:).y]';
        for fk=1:length(f)
            plot(trials(f(fk)).x(:,1),trials(f(fk)).x(:,2)+yoffset*U,'linewidth',.01)
            [val,ind]=findpeaks(abs(trials(f(fk)).f(:,2)),'minpeakheight',5,'npeaks',1);
            time=trials(f(fk)).t;
            try
                fq=find((time>=(time(ind)+onset+fitlower))&(time<=(time(ind)+onset+fitupper)));
            catch
                continue
            end
            quiver(trials(f(fk)).x(fq,1),trials(f(fk)).x(fq,2)+yoffset*U,trials(f(fk)).f(fq,1),trials(f(fk)).f(fq,2),'Color','g')
            plot(trials(f(fk)).x(end,1),trials(f(fk)).x(end,2)+yoffset*U,'rx')
            storeme(U,fk).X=[trials(f(fk)).x(fq,:) trials(f(fk)).v(fq,:) trials(f(fk)).a(fq,:) trials(f(fk)).f(fq,:)];
            Y=gaussianWeightedRegression(catx(:,1),caty,trials(f(fk)).x(fq,1),1000);
            Y(:,1)=Y(:,1)-(Y(1,1)-trials(f(fk)).x(fq(1),2)); %If the "shape" holds, this corrects for starting bias.
            storeme(U,fk).Y=Y;
            [storeme(U,fk).Ts,storeme(U,fk).E,storeme(U,fk).Eb,storeme(U,fk).oT]=cart2model(storeme(U,fk).X,storeme(U,fk).Y);
            storeme(U,fk).trialnum=f(fk)*ones(size(storeme(U,fk).Ts,1),1);
            storeme(U,fk).time=time(fq)-time(fq(1));
            storeme(U,fk).forceson=find((time>=(time(ind)+onset))&(time<=(time(ind)-onset)));
            plot(trials(f(fk)).x(fq,1),storeme(U,fk).Y(:,1)+yoffset*U,'c')
        end
    end
    axis equal

    Eb=vertcat(storeme(:).Eb);
    Eb=[Eb(:,1)+Eb(:,3); Eb(:,2)+Eb(:,4)];
    Ts=vertcat(storeme(:).Ts);
    Ts=[Ts(:,1); Ts(:,2)];
    rT=vertcat(storeme(:).trialnum);
    rT=[rT; rT];
    u=unique(rT); %has the effect of sorting them too

    time=vertcat(storeme(:).time);
    time=[time;time];
    figure(2666)
    clf
    hold on
    for kk=1:length(time)
        plot(Eb(kk),Ts(kk),'.','markersize',1,'color',[.3,time(kk)/.25,0])
    end
    for kk=0:.01:fitupper
        X=Eb(time<=kk);
        Y=Ts(time<=kk);
        W=X\Y;
        minX=min(X);
        maxX=max(X);
        plot([minX,maxX],W*[minX,maxX],'-','color',[0 kk/fitupper 1])
        text(maxX,W*maxX,num2str(kk))
    end

    nT=ceil(log2(length(u)));
    W=zeros(nT+2,3);
    W(1,1)=1;
    W(1,2)=0;
    W(1,3)=mean(abs(Ts-Eb));
    labels{1}='0';

    for U=0:nT
        trainI=rT<=u(min(length(u),ceil(2^U)));
        xtrain=Eb(trainI);
        ytrain=Ts(trainI);
        W(U+2,1)=xtrain\ytrain;
        W(U+2,2)=mean(abs(ytrain-W(U+2,1)*xtrain));
        W(U+2,3)=mean(abs(Ts-W(U+2,1)*Eb));
        labels{U+2}=num2str(min(length(u),ceil(2^U)));
    end
    figure(6)
    clf
    subplot(2,1,1)
    plot(-1:U,W(:,1),'-')
    ylabel('Fit Feedback Gain, unitless')
    set(gca,'xticklabels',[])
    subplot(2,1,2)
    plot(-1:U,W(:,2),'-b',-1:U,W(:,3),'-r')
    ylabel('Mean Error, Nm')
    xlabel('Trials Used to Determine Weight')
    set(gca,'xticklabels',labels)

    baselineCatme=catme;
    save(['../Data/Data_pulse/pulse',num2str(k),'W.mat'],'W','labels','starts','ends','dcats','storeme','baselineCatme')
end



% Specify desired trajectories as replacement of y-component with
% reference via regression on x-component for 50ms following onset

% Convert cartesian space to joint space

% Analytically solve for kp0 and kp1 GAINS though you could simply
% solve... maybe both? Remember that this is essentially a force
% balance: what you should find = k * what you found

% Since these are fast left-divisions: Do trial-by-trial (10 points per
% left division is nothing...) and aggregate plots

