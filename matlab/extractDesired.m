%function extractDesired(name, debug)

name='11'
disp(['Extracting Desired Trajectories for Subject ',name])

%if nargin<2
    debug=0;
%end

load(['../Data/',name,'.mat']);
global kp measuredVals measuredTime x0 fJ getAccel desiredVals desiredTime

params

if 1
    set2dGlobals(params.l1, params.l2, params.origin, params.shoulder,params.mass)
    disp('Mass set')
else
    set2dGlobals(params.l1, params.l2, params.origin, params.shoulder)
end

%Do the extraction on trials where forces were on
lT=length(trials);
for k=1:lT
    length(trials)
    k/lT %#ok<NOPRT>

    inds=trials(k).first:trials(k).last;

    x0=trials(k).x0;
    target=trials(k).target;

    measuredVals=[trials(k).q(inds,:) trials(k).qdot(inds,:) trials(k).qddot(inds,:) trials(k).force(inds,:)];
    measuredTime=trials(k).time(inds)-trials(k).time(trials(k).first);

    if debug
        figure(k)
        clf
        subplot(8,1,1:7)
        hold on
    end

    [T,X]=ode45(@armdynamicsInvertedBurdet,measuredTime,[trials(k).q(inds(1),:)';trials(k).qdot(inds(1),:)']);
    
    desiredTime=T;

    desired.qDesired=X(:,1:2);
    desired.qdotDesired=X(:,3:4);
    desired.xDesired=zeros(length(T),2);
    desired.vDesired=zeros(length(T),2);
    desired.aDesired=zeros(length(T),2);
    desired.time=measuredTime;

    desiredVals=zeros(length(T),8);
    desiredVals(:,1:4)=X;

    for kk=1:length(T)
        desired.xDesired(kk,:)=fkin(desired.qDesired(kk,:));
        desired.vDesired(kk,:)=(fJ(desired.qDesired(kk,:))*desired.qdotDesired(kk,:)')';
    end

    desiredTrajectories(k,1)=desired; %#ok<NASGU>

    [T,X]=extractionReflexHelper(measuredTime,[trials(k).q(inds(1),:)';trials(k).qdot(inds(1),:)'],[ikin(trials(k).des(1,:))' 0 0]');
    desiredTime=T;

    desired.qDesired=X(:,1:2);
    desired.qdotDesired=X(:,3:4);
    desired.xDesired=zeros(length(T),2);
    desired.vDesired=zeros(length(T),2);
    desired.aDesired=zeros(length(T),2);
    desired.time=measuredTime;

    desiredVals=zeros(length(T),8);
    desiredVals(:,1:4)=X;

    for kk=1:length(T)
        desired.xDesired(kk,:)=fkin(desired.qDesired(kk,:));
        desired.vDesired(kk,:)=(fJ(desired.qDesired(kk,:))*desired.qdotDesired(kk,:)')';
    end

    desiredTrajectories(k,2)=desired; %#ok<NASGU>
end

save(['../Data/',name,'extracted.mat'],'desiredTrajectories');