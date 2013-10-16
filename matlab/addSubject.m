%function addSubject(name)
clc
clear all
name='301'

doPlots=1;
%doPlots=0;

newline=sprintf('\n');
%% Get everything loaded and setup
disp(['Loading Data for Subject ',name])
tic

output=load(['../Data/output301.dat']);
input=load(['../Data/input301.dat']);
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
loadT=toc;

%% Smooth everything at once because it's actually faster and more efficient.
tic
s01=size(output,1);
trial=output(:,1);

filtn=128;
filtType='loess';

t=linspace(output(1,13),output(end,13),s01)';
gT=t(2)-t(1);
x=[smooth(output(:,3),filtn,filtType) smooth(output(:,4),filtn,filtType)];
xdiff=x-output(:,3:4);
var(xdiff)

disp('Position data smoothed.')
toc
disp(newline)

v=[gradient(x(:,1)) gradient(x(:,2))]/gT;
a=[gradient(v(:,1)) gradient(v(:,2))]/gT;
disp('State data differentiated.')
toc
disp(newline)

f=[smooth(output(:,9),filtn,filtType) smooth(output(:,10),filtn,filtType)];
fdiff=f-output(:,9:10);

disp('Force data smoothed.')
toc
disp(newline)
smoothT=toc;

%% Our force sensors gradually accrue an offset. Find spots where minimal force is being applied to the handle
tic
calib=find(vmlt(a,.1)&vmlt(v,.005)); %Square roots are unnecessarily slow.
disp('Calibration Points Found.')
toc
disp(newline)

minSpanT=.2; %300 ms
minSpanI=ceil(minSpanT/gT); %Convert to index span

dtimec=[0; diff(t(calib))]; %Time distance since calibration was satisfied
breaks=find(dtimec>.1);

begins=calib(breaks(1:end-1));
ends=calib(breaks(2:end)-1);

c=0;
clear calibclumps
for k=1:length(begins)
    if (ends(k)-begins(k))>=minSpanI
        c=c+1;
        calibclumps(c).inds=begins(k):ends(k);
        calibclumps(c).center=round(mean(calibclumps(c).inds));
    end
end
disp('Calibration Time Segments Determined.')
toc
disp(newline)

if doPlots
    upper=length(t);
    subinds=1:upper;

    figure(1)
    clf
    hold on
    plot(t(subinds),v(subinds,1),'b')
    plot(t(begins(begins<upper)),v(begins(begins<upper),1),'mx')
    plot(t(ends(ends<upper)),v(ends(ends<upper),1),'kx')
    for k=1:length(calibclumps)
        if calibclumps(k).inds(end)<subinds(end)
            plot(t(calibclumps(k).inds),v(calibclumps(k).inds,1),'r.')
        end
    end
end

clumpT=toc;

%% Now deal with rotating the forces to deal with that accrued offset
tic

%Trig is easily the slowest thing we do, store everything we compute
qr=x;
q=t;
sq=q;
cq=q;
frot=f;

for k=1:s01
    qr(k,:)=ikinRobot(x(k,:));
    q=sum(qr(k,:));
    s=sin(q);
    c=cos(q);
    sq(k)=s;
    cq(k)=c;
    frot(k,:)=([c -s;s c]'*frot(k,:)')'; %Rotate forces from roomspace to a consistent robot handle space: q_robot=[0 0]
end
disp('Rotations Performed')
toc
disp(['Rate of rotation: ',num2str((t(k)-t(1))/toc),' secs data/sec computation'])
disp(newline)

trigT=toc;

%% Use interpolation between calibration segments to fix rotated forces
tic

frotfix=frot;

%Deal with data prior to the first clump
mc=mean(f(calibclumps(1).inds,:));
zeroinds=1:calibclumps(1).center-1;
frotfix(zeroinds,:)=ones(length(zeroinds),1)*mc;

for k=1:length(calibclumps)-1
    fitinds=[calibclumps(k).inds calibclumps(k+1).inds];
    px=polyfit(t(fitinds),frot(fitinds,1),1);
    py=polyfit(t(fitinds),frot(fitinds,2),1);

    fixinds=calibclumps(k).center:calibclumps(k+1).center-1;
    frotfix(fixinds,1)=polyval(px,t(fixinds));
    frotfix(fixinds,2)=polyval(py,t(fixinds));
end

%Deal with data after to the last clump
mc=mean(f(calibclumps(end).inds,:));
endinds=calibclumps(end).center:length(t);
frotfix(endinds,:)=ones(length(endinds),1)*mc;

disp('Force Fix Calculated.')
toc
disp(newline)

if doPlots
    figure(2)
    clf
    hold on
    
    upper=100000; %length(t);
    subinds=1:100:upper;
    
    plot(t(subinds),frot(subinds,1),'b')
    fixedRot=frot-frotfix;
    plot(t(subinds),fixedRot(subinds,1),'Color',[.5 .5 .5])
    plot(t(subinds),frotfix(subinds,1),'k')
    for k=1:length(calibclumps)
        if (calibclumps(k).inds(1)>=subinds(1))&&(calibclumps(k).inds(end)<=subinds(end))
            plot(t(calibclumps(k).inds),fixedRot(calibclumps(k).inds,1),'r.')
        end
    end
    plot(t(subinds([1 end])),[0 0],'m-')
    title('Opportunistic Force Recalibration in Robot Handle Space')
    xlabel('Time, s')
    ylabel('Force_X, N')
    legend('Original','Fixed','Difference','Recalibration Clusters')
end

mean(frotfix)
std(frotfix)

forcefixT=toc;

%% Use trig again to rotate the fix into lab-space coordinates
tic

ffix=frotfix;

for k=1:s01
    s=sq(k);
    c=cq(k);
    ffix(k,:)=([c -s;s c]*frotfix(k,:)')'; %Rotate forces from roomspace to a consistent robot handle space: q_robot=[0 0]
end

disp('Force Fix Rotated.')
toc
disp(newline)

if doPlots
    figure(3)
    clf
    hold on
    
    upper=100000 %length(t);
    subinds=1:100:upper;
    
    plot(t(subinds),frot(subinds,1),'b')
    fixed=f-ffix;
    plot(t(subinds),fixed(subinds,1),'Color',[.5 .5 .5])
    plot(t(subinds),ffix(subinds,1),'k')
    for k=1:length(calibclumps)
        if (calibclumps(k).inds(1)>=subinds(1))&&(calibclumps(k).inds(end)<=subinds(end))
            plot(t(calibclumps(k).inds),fixed(calibclumps(k).inds,1),'r.')
        end
    end
    plot(t(subinds([1 end])),[0 0],'m-')
    title('Opportunistic Force Recalibration in Lab Space')
    xlabel('Time, s')
    ylabel('Force_X, N')
    legend('Original','Fixed','Difference','Recalibration Clusters')
end

fixrotT=toc;

%% Set up trial structure and save it

u=unique(trial);

for k=1:length(u)
    inds=find(trial==u(k));
    if length(inds)<(.2/gT) %Ie. Not even 200 ms long...
        continue
    end

    trials(k).rawnum=u(k);
    trials(k).x=x(inds,:);
    trials(k).v=v(inds,:);
    trials(k).a=a(inds,:);
    trials(k).f=fixed(inds,:);
    trials(k).fraw=f(inds,:);
    trials(k).t=t(inds);

    if u(k)>1
        trials(k).origin=input(u(k)-1,[2 3]);
    else
        trials(k).origin=trials(k).x(1,:);
    end
    trials(k).target=input(u(k),[2 3]);
    
    trials(k).early=input(u(k),4);
    trials(k).late=input(u(k),5);
    trials(k).white=input(u(k),6);
    trials(k).updown=sum(input(u(k),4:6));
    trials(k).shape=input(u(k),7);
    trials(k).vision=input(u(k),8);
    
    dir=trials(k).target-trials(k).origin;
    trials(k).direction=atan2(dir(2),dir(1));
    trials(k).dist=norm(dir);
end

total_time=fixrotT+forcefixT+trigT+clumpT+smoothT+loadT

save(['../Data/',name,'.mat'],'trials','params');
%Notice structure array and alignment of data designed for concatenation