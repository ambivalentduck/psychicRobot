%function success=addSubject(name)
clc
clear all
name='500'

disp(['Loading Data for Subject ',name])

output=load(['../Data/output500.dat']);
input=load(['../Data/input_forcetest.dat']);

global fJ getAlpha x0

LEFT=.281;
RIGHT=-.3951;
TOP=.18;
BOTTOM=.68;
origin=[(LEFT+RIGHT)/2,(TOP+BOTTOM)/2];

[l1, l2, shoulder,mass]=getSubjectParams('1');
%Seems reasonable to measure arms, placement. Unreasonable to weigh.
params.l1=l1;
params.l2=l2;
params.shoulder=shoulder;
params.origin=origin;
params.dimensions=2;
params.mass=mass;
set2dGlobals(l1, l2, origin, shoulder,mass)
x0_=x0;

%trial TAB now-zero TAB cursor.X() TAB cursor.Y() TAB velocity.X() TAB velocity.Y() TAB accel.X() TAB accel.Y() TAB force.X() TAB force.Y() TAB sigGain

%f=find(sum(abs(input(1:730,4:6)),2)>0);
f=9:50;
a=f; %unique(input(f,1));

success=zeros(length(a),5);

for k=1:length(a)
    K=a(k)
    fo=find(output(:,1)==K);
    trials(k).rawnum=K;
    trials(k).target=input(K,[2 3]);
    trials(k).origin=input(K-1,[2 3]);
    trials(k).early=input(K,4);
    trials(k).late=input(K,5);
    trials(k).white=input(K,6);
    trials(k).updown=sum(input(K,4:6));
    trials(k).shape=input(K,7);
    trials(k).type=find(input(K,4:6)~=0);
        
    trials(k).time=output(fo,2);
    trials(k).pos=output(fo,[3 4]);
    trials(k).des=output(fo,[11 12]);
    gT=(output(fo(end),2)-output(fo(1),2))/length(fo);
    trials(k).desvel=[gradient(trials(k).des(:,1))./gT gradient(trials(k).des(:,2))./gT];
    if norm(trials(k).pos(1,:)-trials(k).pos(end,:))<.25
        trials(k).long=0;
    else
        trials(k).long=1;
    end

    trials(k).vel=output(fo,[5 6]);
    trials(k).accel=output(fo,[7 8]);
    trials(k).force=output(fo,[9 10]);
    %trials(k).force=[-trials(k).force(:,1) trials(k).force(:,2)];

    trials(k).dist=[0; cumsum(sqrt(sum((trials(k).pos(2:end,:)-trials(k).pos(1:end-1,:)).^2,2)))];
    trials(k).speed=sqrt(sum(trials(k).vel.^2,2));

    trials(k).q=trials(k).pos;
    trials(k).qdot=trials(k).q;
    trials(k).qddot=trials(k).q;
    trials(k).torque=trials(k).force;

    flag=1;
    x0=x0_;
    warning off all

    gT=.01; %gradient(trials(k).time);

    while flag
        for kk=1:length(trials(k).time)
            trials(k).q(kk,:)=ikin(trials(k).pos(kk,:));
            trials(k).qdot(kk,:)=(fJ(trials(k).q(kk,:))\(trials(k).vel(kk,:)'))';
            trials(k).qddot(kk,:)=getAlpha(trials(k).q(kk,:)',trials(k).qdot(kk,:)',trials(k).accel(kk,:)');
            trials(k).torque(kk,:)=((fJ(trials(k).q(kk,:))')*(trials(k).force(kk,:)'))';
        end
        %trials(k).qdot=[gradient(trials(k).q(:,1))./gT gradient(trials(k).q(:,2))./gT];clc
        %trials(k).qddot=[gradient(trials(k).qdot(:,1))./gT gradient(trials(k).qdot(:,2))./gT];
        success(k,1:4)=[k K isreal(trials(k).q) x0(2)];
        if success(k,3)
            flag=0;
        else
            x0(2)=x0(2)-.01; %It's obvious that their joints never go imaginary. Forward shoulder movement is the most likely explaination
        end
    end
    for kk=1:length(trials(k).time)
        if norm(trials(k).pos(kk,:)-fkin(trials(k).q(kk,:))')>.001
            [k kk trials(k).pos(kk,:) trials(k).q(kk,:)]
        end
    end

    warning on all

    trials(k).x0=x0;
    trials(k).targetCat=trials(k).pos(end,1)>0;

    f=find((trials(k).time-trials(k).time(1))<.5,1,'last'); 
    trials(k).first=f;
    trials(k).last=length(trials(k).time);
    
    %Should adjust to LOW velocity/acceleration threshold, not just mean.
    %trials(k).force=trials(k).force-ones(length(trials(k).time),1)*mean(trials(k).force(1:trials(k).first,:));

    success(k,5)=trials(k).first;
end

success %#ok<NOPRT>

save(['../Data/',name,'.mat'],'trials','params');
%Notice structure array and alignment of data designed for concatenation