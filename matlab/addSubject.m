%function success=addSubject(name)
clc
clear all
name='300'

newline=sprintf('\n');
%% Get everything loaded and setup
disp(['Loading Data for Subject ',name])
tic

output=load(['../Data/output300.dat']);
input=load(['../Data/input.dat']);
disp('Data files loaded.')
toc
disp(newline)

global fJ getAlpha x0

LEFT=.281;
RIGHT=-.3951;
TOP=.18;
BOTTOM=.68;
origin=[(LEFT+RIGHT)/2,(TOP+BOTTOM)/2];

[l1, l2, shoulder,mass]=getSubjectParams(name);
%Seems reasonable to measure arms, placement. Unreasonable to weigh.
params.l1=l1;
params.l2=l2;
params.shoulder=shoulder;
params.origin=origin;
params.dimensions=2;
params.mass=mass*.4535; %lbs to kg
set2dGlobals(l1, l2, origin, shoulder,mass)
x0_=x0;
toc
disp(newline)

%% Smooth everything at once because it's actually faster and more efficient.
tic
s01=size(output,1);
trial=output(:,1);

filtn=128;
filtType='loess';

t=linspace(output(1,13),output(end,13),s01)';
gT=t(2)-t(1);
x=[smooth(output(:,3),filtn,filtType) smooth(output(:,4),filtn,filtType)];
disp('Position data smoothed.')
toc
disp(newline)

v=[gradient(x(:,1)) gradient(x(:,2))]/gT;
a=[gradient(v(:,1)) gradient(v(:,2))]/gT;
disp('State data differentiated.')
toc
disp(newline)

f=[smooth(output(:,9),filtn,filtType) smooth(output(:,10),filtn,filtType)];
disp('Force data smoothed.')
toc
disp(newline)

%% Our force sensors gradually accrue an offset. Find spots where minimal force is being applied to the handle
tic
calib=find(vmlt(a,.1)&vmlt(v,.005)); %Square roots are unnecessarily slow.
disp('Calibration Points Found.')
toc
ltoc=toc;
disp(newline)

minSpanT=.2; %300 ms
minSpanI=ceil(minSpanT/gT); %Convert to index span

dtimec=[0; diff(t(calib))]; %Time distance since calibration was satisfied
breaks=find(dtimec>.1)

begins=calib(breaks(1:end-1));
ends=calib(breaks(2:end)-1);

c=0;
clear calibclumps
for k=1:length(begins)
    if (ends(k)-begins(k))>=minSpanI
        c=c+1;
        calibclumps(c).inds=begins(k):ends(k);
    end
end

upper=100000;
subinds=1:upper;

figure(1)
clf
hold on
plot(t(subinds),v(subinds,1),'b')
plot(t(begins(begins<upper)),v(begins(begins<upper),1),'mx')
plot(t(ends(ends<upper)),v(ends(ends<upper),1),'kx')
for k=1:length(lumps)
    if lumps(k).inds(end)<subinds(end)
        plot(t(calibclumps(k).inds),v(calibclumps(k).inds,1),'r.')
    end
end


% oops=find((~dtrialc)&(dtimec>.1)) %Places where our calibration point finder additional "pauses" in a trial.
% 
% figure(1)
% clf
% hold on
% N=20;
% ti=find(trial==trial(calib(oops(N))));
% o=calib(oops(trial(calib(oops))==trial(calib(oops(N)))));
% plot(t(ti),v(ti),'b')
% ti2=((vmlt(a,.1)&vmlt(v,.005))&(trial==trial(calib(oops(N)))));
% plot(t(ti2),v(ti2),'r.')
% plot(t(o),v(o),'kx')


%% Now deal with rotating the forces to deal with that accrued offset

return 
qr=x;
q=t;
sq=q;
cq=q;
frot=f;
frot2=frot;

for k=1:s01
    qr(k,:)=ikinRobot(x(k,:)); %Trig is easily the slowest thing we do, store everything we compute
    q=sum(qr(k,:));
    s=sin(q);
    c=cos(q);
    sq(k)=s;
    cq(k)=c;
    frot(k,:)=([c -s;s c]'*frot(k,:)')'; %Rotate forces from roomspace to a consistent robot handle space: q_robot=[0 0]
end
disp('Rotations Performed')
toc
ltoc=toc-ltoc;
disp(['Rate of rotation: ',num2str((t(k)-t(1))/ltoc),' secs data/sec computation'])
disp(newline)


%% Later stuff

a=find(sum(abs(input(:,4:6)),2)>0);
%a=1:66;

success=zeros(length(a),5);

for k=1:length(a)
    K=a(k)
    fo=find(output(:,1)==K);
    trials(k).rawnum=K;
    trials(k).target=input(K,[2 3]);
    if K>1
        trials(k).origin=input(K-1,[2 3]);
    else
        trials(k).origin=input(fo(1),[3 4]);
    end

    trials(k).early=input(K,4);
    trials(k).late=input(K,5);
    trials(k).white=input(K,6);
    trials(k).updown=sum(input(K,4:6));
    trials(k).shape=input(K,7);
    trials(k).type=find(input(K,4:6)~=0);
    trials(k).vision=input(K,8);

    %trials(k).time=output(fo,13);
    trials(k).time=linspace(output(fo(1),13),output(fo(end),13),length(fo));
    gT=mean(gradient(trials(k).time));
    
    trials(k).pos=output(fo,[3 4]);
    trials(k).des=output(fo,[11 12]);
    
    trials(k).desvel=[gradient(trials(k).des(:,1))./gT gradient(trials(k).des(:,2))./gT];
    if norm(trials(k).pos(1,:)-trials(k).pos(end,:))<.25
        trials(k).long=0;
    else
        trials(k).long=1;
    end


    gT=mean(gradient(trials(k).time));
    trials(k).vel=[gradient(trials(k).pos(:,1))./gT gradient(trials(k).pos(:,2))./gT];
    trials(k).accel=[gradient(trials(k).vel(:,1))./gT gradient(trials(k).vel(:,2))./gT];
    %trials(k).vel=output(fo,[5 6]);
    %trials(k).accel=output(fo,[7 8]);
    trials(k).force=output(fo,[9 10]);
    
    for kk=1:length(trials(k).time)
        trials(k).qRobot(kk,:)=ikinRobot(trials(k).pos(kk,:));
        trials(k).qdotRobot(kk,:)=(robotfJ(trials(k).qRobot(kk,:))\(trials(k).vel(kk,:)'))';
        trials(k).qddotRobot(kk,:)=robotAlpha(trials(k).qRobot(kk,:)',trials(k).qdotRobot(kk,:)',trials(k).accel(kk,:)');
        trials(k).torqueRobot(kk,:)=((robotfJ(trials(k).qRobot(kk,:))')*(trials(k).force(kk,:)'))';
    end
        
    trials(k).dist=[0; cumsum(sqrt(sum((trials(k).pos(2:end,:)-trials(k).pos(1:end-1,:)).^2,2)))];
    trials(k).speed=sqrt(sum(trials(k).vel.^2,2));

    trials(k).q=trials(k).pos;
    trials(k).qdot=trials(k).q;
    trials(k).qddot=trials(k).q;
    trials(k).torque=trials(k).force;

    flag=1;
    x0=x0_;
    warning off all

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