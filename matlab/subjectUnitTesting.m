%function success=addSubject(name)
clc
%clear all
name='300'

disp(['Loading Data for Subject ',name])

%output=load(['../Data/output',name,'.dat']);
input=load(['../Data/input.dat']);

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

global kpgain
kpgain=1;

%% We've had issues with irregular time traces
% For two undisturbed movements and a kicked movement, plot velocity traces
% from various methods.
figure(1)
clf
MOVES=[8 600 find(input(:,5)>0,1,'first')];
lM=length(MOVES);
for k=1:lM

    f=find(output(:,1)==MOVES(k));
    x=output(f,[3 4]);
    subplot(4,lM,k)
    plot(x(:,1),x(:,2),'.')
    axis equal
    title(['Movement ',num2str(MOVES(k))])
    if(k==1)
        ylabel('Position Trace')
    end
    xlabel('Robot Space, m')

    subplot(4,lM,k+lM)
    t=output(f,13);
    gT=gradient(t);
    u=unique(gT)
    h=hist(gT,u);
    bar(1000*u,h)
    if(k==1)
        ylabel('Histogram: Clock Skips?')
    end
    xlabel('Change in reported time, ms')

    subplot(2,lM,k+lM)
    hold on
    plot(t,output(f,5),'g.')

    f=find(gT~=0);
    gT=gT(f);
    tf=t(f);
    xf=x(f,:);
    vf=[gradient(xf(:,1))./gT gradient(xf(:,2))./gT];
    plot(tf,vf(:,1),'bo')
    tl=linspace(t(1),t(end),length(t));
    vl=[gradient(x(:,1)) gradient(x(:,2))]/((t(end)-t(1))/(length(t)-1));
    plot(t,vl(:,1),'r.')

    if(k==1)
        ylabel('Velocity Traces')
    end
    if(k==lM)
        legend('Felix''s Value','Throw out grad(t)==0','linearly spaced t')
    end
    xlabel('Time, s')

end
suplabel('Before filtering','t')

%% The velocity traces above are a disaster, so smooth. But how much?
smoothMeth='loess';
figure(2)
clf

N=0:11;
N=2.^N;
Y=zeros(length(t),length(N));
V=Y;
A=V;
gT=mean(diff(tl));
for n=1:length(N)
    Y(:,n)=smooth(x(:,1),N(n),smoothMeth);
    V(:,n)=gradient(Y(:,n))/gT;
    A(:,n)=gradient(V(:,n))/gT;
end
SPACER=20;
[p,q]=meshgrid(log2(N),t(1:SPACER:end)-t(1));

subplot(2,1,1)
mesh(p,q,Y(1:SPACER:end,:))
xlabel('Filter length, log2')
ylabel('Time, s')
zlabel('Robot X, m')

subplot(2,1,2)
mesh(p,q,A(1:SPACER:end,:))
xlabel('Filter length, log2')
ylabel('Time, s')
zlabel('Robot vel_x, m/s')

%% Using a filter length from what we see above, redo fig 1
figure(3)
clf
MOVES=[8 600 find(input(:,5)>0,1,'first')];
lM=length(MOVES);
for k=1:lM

    f=find(output(:,1)==MOVES(k));
    xraw=output(f,[3 4]);
    filtn=128;
    filttype='loess';
    x=[smooth(xraw(:,1),filtn,filttype) smooth(xraw(:,2),filtn,filttype)];
    subplot(4,lM,k)
    plot(x(:,1),x(:,2),'.')
    axis equal
    title(['Movement ',num2str(MOVES(k))])
    if(k==1)
        ylabel('Position Trace')
    end
    xlabel('Robot Space, m')

    subplot(4,lM,k+lM)
    t=output(f,13);
    gT=gradient(t);
    u=unique(gT)
    h=hist(gT,u);
    bar(1000*u,h)
    if(k==1)
        ylabel('Histogram: Clock Skips?')
    end
    xlabel('Change in reported time, ms')

    subplot(2,lM,k+lM)
    hold on
    plot(t,output(f,5),'g.')

    f=find(gT~=0);
    gT=gT(f);
    tf=t(f);
    xf=x(f,:);
    vf=[gradient(xf(:,1))./gT gradient(xf(:,2))./gT];
    plot(tf,vf(:,1),'bo')
    tl=linspace(t(1),t(end),length(t));
    vl=[gradient(x(:,1)) gradient(x(:,2))]/((t(end)-t(1))/(length(t)-1));
    plot(t,vl(:,1),'r.')

    if(k==1)
        ylabel('Velocity Traces')
    end
    if(k==lM)
        legend('Felix''s Value','Throw out grad(t)==0','linearly spaced t')
    end
end

%% Now that the velocity traces aren't awful: xvaf
figure(4)
clf
filtn=128;
filttype='loess';
pad=filtn; %remove edge effects
SPACER=15;

MOVES=[8 600 find(input(:,5)>0,1,'first')];
lM=length(MOVES);
for k=1:lM
    f=find(output(:,1)==MOVES(k));
    fpad=f(1)-pad:f(end)+pad;
    xraw=output(fpad,[3 4]);
    t=output(f,13);
    tl=linspace(t(1),t(end),length(t));
    mgtl=tl(2)-tl(1);
    x=[smooth(xraw(:,1),filtn,filttype) smooth(xraw(:,2),filtn,filttype)];
    x=x(pad+1:end-pad,:);

    subplot(4,lM,k)
    hold on
    plot(x(1:SPACER:end,1),x(1:SPACER:end,2),'.',xraw(1:SPACER:end,1),xraw(1:SPACER:end,2),'o')
    axis equal
    title(['Movement ',num2str(MOVES(k))])
    if(k==1)
        ylabel('Position Trace')
    end
    xlabel('Robot Space, m')

    subplot(4,lM,k+lM)
    v=[gradient(x(:,1)) gradient(x(:,2))]/mgtl;
    vraw=output(f,[5 6]);
    plot(tl(1:SPACER:end),v(1:SPACER:end,1),'.',tl(1:SPACER:end),vraw(1:SPACER:end,1),'o')
    if(k==1)
        ylabel('Velocity_x Trace, m/s')
    end

    subplot(4,lM,k+2*lM)
    a=[gradient(v(:,1)) gradient(v(:,2))]/mgtl;
    araw=output(f,[7 8]);
    plot(tl(1:SPACER:end),a(1:SPACER:end,1),'.',tl(1:SPACER:end),araw(1:SPACER:end,1),'o')
    if(k==1)
        ylabel('Accel_x Trace, m/s^2')
    end

    subplot(4,lM,k+3*lM)
    fraw=output(fpad,[9 10]);
    f=[smooth(fraw(:,1),filtn,filttype) smooth(fraw(:,2),filtn,filttype)];
    f=f(pad+1:end-pad,:);
    fraw=fraw(pad+1:end-pad,:);

    plot(tl(1:SPACER:end),f(1:SPACER:end,1),'.',tl(1:SPACER:end),fraw(1:SPACER:end,1),'o')
    if(k==1)
        ylabel('Force_x Trace, N')
    end

    subplot(4,lM,k)

    kpgain=1;
    y=extract(tl,[x v a f],'reflex');
    plot(y(:,1),y(:,2),'r')
    if k==lM
        legend('Filtered','Raw','Extracted')
    end
end
suplabel('Filtering Outcome','t')

%% The extraction is even pretty reasonable, BUT force offset remains.
figure(5)
clf

figure(6)
clf

filtn=128;
filttype='loess';
pad=filtn*3; %remove edge effects
SPACER=20;
qscale=.5;

MOVES=[8 600 find(input(:,5)>0,1,'first')];
lM=length(MOVES);
for k=1:lM
    figure(5)
    
    f=find(output(:,1)==MOVES(k));
    indices=f;
    fpad=f(1)-pad:f(end)+pad;
    xraw=output(fpad,[3 4]);
    t=output(f,13);
    t0=t(1);
    t=t-t0;
    tl=linspace(t(1),t(end),length(t));
    mgtl=tl(2)-tl(1);
    tpad=-pad*mgtl:mgtl:tl(end)+pad*mgtl;
    xrs=[smooth(xraw(:,1),filtn,filttype) smooth(xraw(:,2),filtn,filttype)];
    vpad=[gradient(xrs(:,1)) gradient(xrs(:,2))]/mgtl;
    apad=[gradient(vpad(:,1)) gradient(vpad(:,2))]/mgtl;
    x=xrs(pad+1:end-pad,:);
    v=[gradient(x(:,1)) gradient(x(:,2))]/mgtl;
    a=[gradient(v(:,1)) gradient(v(:,2))]/mgtl;
    fraw=output(fpad,[9 10]);
    f=[smooth(fraw(:,1),filtn,filttype) smooth(fraw(:,2),filtn,filttype)];
    frot=f;
    f=f(pad+1:end-pad,:);
    
    qr=xrs;
    for n=1:length(fpad)
        qr(n,:)=ikinRobot(xrs(n,:));
        q=sum(qr(n,:));
        s=sin(q);
        c=cos(q);
        frot(n,:)=([c -s;s c]*frot(n,:)')';
    end
    frup=frot(pad+1:end-pad,:);
    calib=find((vecmag(apad)<.3)&(vecmag(vpad)<.005));
    
    subplot(3,lM,k)
    hold on
    quiver(x(1:SPACER:end,1),x(1:SPACER:end,2),f(1:SPACER:end,1),f(1:SPACER:end,2),qscale,'Color',[.5 .5 .5])
    plot(x(1:SPACER:end,1),x(1:SPACER:end,2),'b')
    plot(xrs(calib,1),xrs(calib,2),'rx')
    if(k==1)
        ylabel('Position trace with forces')
    end
    
    subplot(3,lM,k+lM)
    hold on
    quiver(tl(1:SPACER:end)',x(1:SPACER:end,1),f(1:SPACER:end,1),f(1:SPACER:end,2),qscale,'Color',[.5 .5 .5])
    quiver(tpad(1:SPACER:end)',xrs(1:SPACER:end,1),frot(1:SPACER:end,1),frot(1:SPACER:end,2),qscale,'c')
    plot(tl(1:SPACER:end),x(1:SPACER:end,1),'b')
    plot(tpad(calib)',xrs(calib,1),'r.')
    xlabel('Time, s')
    if(k==1)
        ylabel('Forces on pos_x Trace, m')
    end
    
    subplot(3,lM,k+2*lM)
    hold on
    px=polyfit(xrs(calib,1),frot(calib,1),1);
    py=polyfit(xrs(calib,1),frot(calib,2),1);
    
    frot2=frot;
    frot2(:,1)=frot(:,1)-polyval(px,xrs(:,1));
    frot2(:,2)=frot(:,2)-polyval(py,xrs(:,1));

    fcomp=frot;
    for n=1:length(fpad)
        q=sum(qr(n,:));
        s=sin(q);
        c=cos(q);
        fcomp(n,:)=([c -s;s c]'*frot2(n,:)')';
    end
    fcompup=fcomp(pad+1:end-pad,:);
    
    quiver(tpad(1:SPACER:end)',xrs(1:SPACER:end,1),fcomp(1:SPACER:end,1),fcomp(1:SPACER:end,2),qscale,'Color',[.5 .5 .5])
    quiver(tpad(1:SPACER:end)',xrs(1:SPACER:end,1),frot2(1:SPACER:end,1),frot2(1:SPACER:end,2),qscale,'c')
    plot(tl(1:SPACER:end),x(1:SPACER:end,1),'b')
    plot(tpad(calib)',xrs(calib,1),'r.')
    xlabel('Time, s')
    if k==lM
        legend('Smoothed','In Robot Coordinates','Movement Trace','Subset for Force Calibration')
    end
    if(k==1)
        ylabel('Corrected forces on pos_x Trace, m')
    end
    
    figure(6)
    
    yraw=extract(tl,output(indices,3:10),'reflex');
    ysmooth=extract(tl,[x v a f],params,'reflex');
    ysmoothfcorrect=extract(tl,[x v a fcompup],params,'reflex');
    
    subplot(2,lM,k)
    hold on
    plot(x(1:SPACER:end,1),x(1:SPACER:end,2),'b')
    quiver(x(1:SPACER:end,1),x(1:SPACER:end,2),f(1:SPACER:end,1),f(1:SPACER:end,2),qscale,'Color',[.5 .5 .5])
    plot(yraw(:,1),yraw(:,2),'r')
    plot(ysmooth(:,1),ysmooth(:,2),'c')
    plot(ysmoothfcorrect(:,1),ysmoothfcorrect(:,2),'k')
    axis equal
    if k==1
        ylabel('Position Traces')
    end
    subplot(2,lM,k+lM)
    hold on
    plot(tl,vecmag(v),'b')
    plot(tl,vecmag(yraw(:,3:4)),'r')
    plot(tl,vecmag(ysmooth(:,3:4)),'c')
    plot(tl,vecmag(ysmoothfcorrect(:,3:4)),'k')
    
    if k==lM
        legend('Smoothed Recorded','Extracted no smooth','Extracted smooth','Extracted smoothed forces corrected')
    end
    if k==1
        ylabel('Speed Traces')
    end
    xlabel('Time, s')


    
end
figure(5)
suplabel('Fixing the Force Offset','t')


