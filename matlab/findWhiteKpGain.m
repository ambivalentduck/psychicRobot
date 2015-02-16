function findWhiteKpGain(k)

load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])

GWR=2000;

%% Categorize by start/end pair
figure(1)
clf
hold on
clustermeS=zeros(length(trials)-1,1); %#ok<NODEF>
clustermeE=clustermeS;
for t=1:length(trials)-1
    plot(trials(t+1).x([1 end],1),trials(t+1).x([1 end],2),'.')
    clustermeS(t)=trials(t+1).x(1,1);
    clustermeE(t)=trials(t+1).x(end,1);
end
[trash,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
means=sort(means);
plot(means,.5*[1 1 1],'rx')
axis equal

starts=zeros(length(trials),1);
ends=starts;
for t=1:length(trials)
    [trash,starts(t)]=min(abs(trials(t).x(1,1)-means));
    [trash,ends(t)]=min(abs(trials(t).x(end,1)-means));
    trialInfo(t).startcat=starts(t);
    trialInfo(t).endcat=ends(t);
end
starts(1)=-inf; %never use this reach...
trialInfo(1).startcat=-inf;

reachcat=10*starts+ends;
urc=unique(reachcat);
urc=urc(2:end);

%% Use undisturbed examples for a reference trajectory for each start/end pair.
dcats=[trials.disturbcat];

nclean=3;
clean=0*dcats;
for kk=1:nclean
    clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
end
clean=(clean==0);

for t=1:length(trials)
    trialInfo(t).clean=clean(t);
    trialInfo(t).reachcat=find(reachcat(t)==urc);
end

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
    Y=gaussianWeightedRegression(catx(1,:)',caty',X,GWR);
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
    Y=gaussianWeightedRegression(catx(1,:)',caty',X,GWR);
    plot(catx(1,:),caty(2,:)+yoffset*U,'b.','markersize',.01)
    plot(X,Y(:,2)+yoffset*U,'g')
end

%% Shift force trace relative to position trace due to sensor lead.
load forcelag.mat
lag=round(forcelag(k)/5);
f_all=vertcat(trials.f);

for T=3:length(trials)-2 %Quietly fail to fix the first two and last two reaches because they don't matter and are undisturbed
    trials(T).f=f_all(trials(T).rawinds+lag,:); %#ok<*AGROW>
end

%% We'll need to know the onset of pulse forces later
figure(5)
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

figure(6)
clf
hold on

fitlower=.0;
fitupper=.5;
for U=1:length(urc)
    f=find((dcats>0)&(dcats<5)&(reachcat'==(urc(U))));
    catx=[catme(U,:).x]';
    caty=[catme(U,:).y]';
    for fk=1:length(f)
        plot(trials(f(fk)).x(:,1),trials(f(fk)).x(:,2)+yoffset*U,'linewidth',.01)
        [val,ind]=findpeaks(abs(trials(f(fk)).f(:,2)),'minpeakheight',5,'npeaks',1);

        time=trials(f(fk)).t;
        trialInfo(f(fk)).kPeak=ind;
        trialInfo(f(fk)).forceinds=find((time>=(time(ind)+onset))&(time<=(time(ind)-onset)));
    end
end


%% Use a helper function to [x v a], [xd vd ad], and F into mass and viscoelastic ratios
set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

%This code is super straight-forward: error and error velocity vs typical
%error and typical error velocity. F=k1(Da+C)+k2(Ze)
% inertia and coriolis are linear in body mass
% k2 is essentially a stiffness adjustment

fitlower=.2/.005;
fitupper=.4/.005;
indspan=fitlower:fitupper;
for U=1:length(urc)
    f=find((dcats==5)&(reachcat'==(urc(U))));
    catx=[catme(U,:).x]';
    caty=[catme(U,:).y]';
    for fk=1:length(f)
        kk=f(fk);

        onset=find(vecmag(trials(kk).v)>.05,1,'first');
        start=max(onset-35,1);
        time=trials(kk).t(start+indspan);
        time=time-time(1);

        %Find typical reach-perpendicular p,v,a.
        Y=gaussianWeightedRegression(catx(:,1),caty,trials(kk).x(start+indspan,1),GWR);
        %Y(:,1)=Y(:,1)-(Y(1,1)-trials(f(fk)).x(start+indspan(1),2)); %If the "shape" holds, this corrects for starting bias.
        
        X=[trials(kk).x(start+indspan,:) trials(kk).v(start+indspan,:) trials(kk).a(start+indspan,:) trials(kk).f(start+indspan,:)];

        [storeme(U,fk).TMP,storeme(U,fk).TSP,storeme(U,fk).TNP]=cart2model(X,Y);
        storeme(U,fk).time=time;
        storeme(U,fk).X=X;
        storeme(U,fk).Y=Y;
    end
end
axis equal

TMP=vertcat(storeme(:).TMP);
TSP=vertcat(storeme(:).TSP);
TNP=vertcat(storeme(:).TNP);
fit=[TMP(:) TSP(:)]\TNP(:)
time=vertcat(storeme(:).time);
indtime=floor(time*200); %indices

figure(5000)
clf
hold on
X=vertcat(storeme(:).X);
Y=vertcat(storeme(:).Y);
plot(X(:,6),X(:,8),'b.')
plot(Y(:,3),X(:,8),'g.')
plot(X(:,5),X(:,7),'r.')

figure(289)
clf
for TN=1:2
    for T=0:(fitupper-fitlower)
        f=find(indtime==T);
        if isempty(f)
            continue
        end
        ft=[TMP(f,TN), TSP(f,TN)]\TNP(f,TN);
        subplot(2,2,TN)
        hold on
        plot(T*5,ft(1),'r.')
        plot(T*5,ft(2),'b.')
        subplot(2,2,2+TN)
        hold on
        err=TNP(f,TN)-[TMP(f,TN) TSP(f,TN)]*ft;
        plot(T*5,err,'k.')
    end
end

error


massgain=fit(1);
kpgain=fit(2);

save(['../Data/Data_pulse/pulse',num2str(k),'W.mat'],'massgain','kpgain','trialInfo','means','trials','params')
