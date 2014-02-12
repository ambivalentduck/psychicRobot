clc
clear all

global kp0gain kp1gain

for k=4 %:4
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

    if 3 %k==2
        warning off all
        %No clue why, but the input file for #3 is whack. Just infer it.
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
    end

    % Use unperturbed examples for a reference trajectory for each start/end pair.
    dcats=[trials.disturbcat];

    figure(500)
    clf
    plot(dcats0,dcats+.3*rand(size(dcats)),'.')

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
            if krp<=10
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
    fitlower=.07;
    fitupper=.17;
    for U=1:length(urc)
        f=find((dcats>0)&(dcats<5)&(reachcat'==(urc(U))));
        catx=[catme(U,:).x]';
        caty=[catme(U,:).y]';
        for fk=1:length(f)
            plot(trials(f(fk)).x(:,1),trials(f(fk)).x(:,2)+yoffset*U,'linewidth',.01)
            [val,ind]=findpeaks(abs(trials(f(fk)).f(:,2)),'minpeakheight',5,'npeaks',1);
            time=trials(f(fk)).t;
            fq=find((time>=(time(ind)+onset+fitlower))&(time<=(time(ind)+onset+fitupper)));
            quiver(trials(f(fk)).x(fq,1),trials(f(fk)).x(fq,2)+yoffset*U,trials(f(fk)).f(fq,1),trials(f(fk)).f(fq,2),'Color','g')
            plot(trials(f(fk)).x(end,1),trials(f(fk)).x(end,2)+yoffset*U,'rx')
            storeme(U,fk).X=[trials(f(fk)).x(fq,:) trials(f(fk)).v(fq,:) trials(f(fk)).a(fq,:) trials(f(fk)).f(fq,:)];
            Y=gaussianWeightedRegression(catx(:,1),caty,trials(f(fk)).x(fq,1),1000);
            Y(:,1)=Y(:,1)-(Y(1,1)-trials(f(fk)).x(fq(1),2)); %If the "shape" holds, this corrects for starting bias.
            storeme(U,fk).Y=Y;
            [storeme(U,fk).Ts,storeme(U,fk).E,storeme(U,fk).Eb,storeme(U,fk).oT]=cart2model(storeme(U,fk).X,storeme(U,fk).Y);
            plot(trials(f(fk)).x(fq,1),storeme(U,fk).Y(:,1)+yoffset*U,'c')
        end
    end
    axis equal

    figure(6)
    clf
    hold on
    colors='rgbcmy';
    for U=1:length(urc)
        for kk=1:size(storeme,2)
            if isempty(storeme(U,kk).Ts)
                continue
            end
            BOUT=storeme(U,kk).Eb;
            SNE=storeme(U,kk).Ts;
            Gb=[BOUT(:,[1 3]);BOUT(:,[2 4])]\[SNE(:,1); SNE(:,2)];
            %Kx=[storeme(U,kk).K0(:,1) storeme(U,kk).K1(:,1); storeme(U,kk).K0(:,2) storeme(U,kk).K1(:,2)];
            %Sy=[storeme(U,kk).S(:,1);storeme(U,kk).S(:,2)];
            %weights=Kx\Sy;
            plot(Gb(1),Gb(2),[colors(U),'.'])
            %storeme(U,kk).Kx=Kx';
            %storeme(U,kk).Sy=Sy';
            %storeme(U,kk).W=weights;
        end
        E=vertcat(storeme(U,:).E);
        TR=vertcat(storeme(U,:).Ts);
        Ta=vertcat(storeme(U,:).oT);
        [Ta(:,[1 2 5 6]) E]\Ta(:,4)
        Ksm=(E(:,1:4)\TR)'
        continue
        Kx=[storeme(U,:).Kx]';
        Sy=[storeme(U,:).Sy]';
        weights=Kx\Sy;
        plot(weights(1),weights(2),[colors(U),'x'],'markersize',10)
    end
    w=[storeme.W];
    m=[w(1,:)' 0*w(1,:)'+1]\w(2,:)'
    forrange=[storeme.W]';
    wx=min(forrange(:,1)):.1:max(forrange(:,1))
    plot(wx,m(1)*wx+m(2),'k-')
    xlabel('Kp0 Gain')
    ylabel('Kp1 Gain')
    title('Fit to first 100ms of pulse. Colors are direction/length pairs')

    figure(7)
    clf
    hold on
    yoffset=.3;
    xoffset=.4;
    for U=1:length(urc)
        f=find((dcats>0)&(dcats<5)&(reachcat'==(urc(U))));
        catx=[catme(U,:).x]';
        caty=[catme(U,:).y]';
        Kx=[storeme(U,:).Kx]';
        Sy=[storeme(U,:).Sy]';
        weights=Kx\Sy;
        kp0gain=max(.1,weights(1)); %Lower bound because...srsly?
        kp1gain=max(.1,weights(2));

        for fk=1:length(f)
            plot(trials(f(fk)).x(:,1)+xoffset*fk,trials(f(fk)).x(:,2)+yoffset*U,'linewidth',.01)
            plot(trials(f(fk)).x(end,1)+xoffset*fk,trials(f(fk)).x(end,2)+yoffset*U,'rx')
            X=[storeme(U,fk).K0(:,1) storeme(U,fk).K1(:,1); storeme(U,fk).K0(:,2) storeme(U,fk).K1(:,2)];
            Y=[storeme(U,fk).S(:,1);storeme(U,fk).S(:,2)];
            plot(storeme(U,fk).X(:,1)+xoffset*fk,storeme(U,fk).Y(:,1)+yoffset*U,'g')
            weights=X\Y;
            if fk<=inf %5
                %[val,ind]=findpeaks(abs(trials(f(fk)).f(:,2)),'minpeakheight',5,'npeaks',1);
                onset=find(vecmag(trials(f(fk)).v)>.1,1,'first');
                start=max(onset-10,1);
                xvaf=[trials(f(fk)).x trials(f(fk)).v trials(f(fk)).a trials(f(fk)).f];
                y=extract(trials(f(fk)).t(start:end)',xvaf(start:end,:),'reflex');
                trials(f(fk)).y=y;
                trials(f(fk)).ygains=[kp0gain kp1gain];
                plot(y(:,1)+xoffset*fk,y(:,2)+yoffset*U,'r')
                drawnow
            end
        end
    end
    axis equal
    legend('Hand','End of Reach','Presumed Intended for Fit','Extracted Intended')
    title('Testing Fit Quality')

    save(['../Data/Data_pulse/pulse',num2str(k),'y.mat'],'params','trials','storeme')
end



% Specify desired trajectories as replacement of y-component with
% reference via regression on x-component for 50ms following onset

% Convert cartesian space to joint space

% Analytically solve for kp0 and kp1 GAINS though you could simply
% solve... maybe both? Remember that this is essentially a force
% balance: what you should find = k * what you found

% Since these are fast left-divisions: Do trial-by-trial (10 points per
% left division is nothing...) and aggregate plots

