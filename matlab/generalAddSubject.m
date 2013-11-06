function trials=generalAddSubject(name,traw,xvafraw,trial,params)

doPlots=1;

newline=sprintf('\n');

%% Get everything loaded and setup
disp(['Loading Data for Subject ',name])
tic

set2dGlobals(params.l1, params.l2, params.origin, params.shoulder,params.mass) %Make sure mass actually in kg

toc
disp(newline)
loadT=toc;

%% Smooth everything at once because it's actually faster and more efficient.
tic
s01=length(traw);

filtn=128;
filtType='loess';

t=linspace(traw(1),traw(end),s01)';
gT=t(2)-t(1);
x=[smooth(xvafraw(:,1),filtn,filtType) smooth(xvafraw(:,2),filtn,filtType)];
xdiff=x-xvafraw(:,1:2);
var(xdiff)

disp('Position data smoothed.')
toc
disp(newline)

v=[gradient(x(:,1)) gradient(x(:,2))]/gT;
a=[gradient(v(:,1)) gradient(v(:,2))]/gT;
disp('State data differentiated.')
toc
disp(newline)

f=[smooth(xvafraw(:,7),filtn,filtType) smooth(xvafraw(:,8),filtn,filtType)];
fdiff=f-xvafraw(:,7:8);
var(fdiff)

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

c=1;
calibclumps(1).inds=1;
calibclumps(1).center=1;

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
frot=[f(:,1) f(:,2)];

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

%Show key issue: force accrual in time has a steady slope due to force
%sensor error

filtn=60/gT; %1 minute
f2=[smooth(frot(:,1),filtn,filtType) smooth(frot(:,2),filtn,filtType)];
frot=frot-f2;
frot_detrend=frot;

disp('Force Fix Calculated.')
toc
disp(newline)

forcefixT=toc;

%% Use trig again to rotate the fix into lab-space coordinates
tic

f_detrend=f;

for k=1:s01
    s=sq(k);
    c=cq(k);
    f_detrend(k,:)=([c -s;s c]*frot_detrend(k,:)')';
end

fixed=f_detrend; %-ffix; Cluster analysis is really failing here, but who cares...

disp('Force Fix Rotated.')
toc
disp(newline)

if doPlots
    figure(3)
    clf
    hold on
    
    upper=100000; %length(t);
    subinds=1:100:min(upper,length(t));
    
    plot(t(subinds),f(subinds,1),'b')
    plot(t(subinds),fixed(subinds,1),'Color',[.5 .5 .5])
    ffix=f-f_detrend;
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
    trials(k).rawinds=inds;
    trials(k).rawnum=u(k);
    trials(k).x=x(inds,:);
    trials(k).v=v(inds,:);
    trials(k).a=a(inds,:);
    trials(k).f=fixed(inds,:);
    trials(k).fraw=f(inds,:);
    trials(k).t=t(inds);
end

total_time=fixrotT+forcefixT+trigT+clumpT+smoothT+loadT %#ok<NOPRT,NASGU>