%function y=extractLumps(t,xvaf,params,armdynamics)

close all
clear all

name='300'
load(['../Data/',name,'.mat'])

global measuredVals measuredTime errorVals errorTime fJ getAlpha kpgain

%% Setup our globals and measured vals

set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

k=30;
t=trials(k).t';
xvaf=[trials(k).x trials(k).v trials(k).a trials(k).f];
kpgain=1;

measuredVals=xvaf;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)');
    torque=-fJq'*xvaf(k,7:8)';
    measuredVals(k,:)=[q' qdot' qddot' torque'];
end
q0=measuredVals(1,1:4);

calib=find(vmlt(xvaf(:,3:4),.027)); %Square roots are unnecessarily slow.
dtimec=[0 diff(t(calib))]; %Time distance since calibration was satisfied
breaks=find(dtimec>.1);
if isempty(calib)
    start=200;
else
    if isempty(breaks)
        start=calib(end)+1;
    else
        start=calib(breaks(1)-1);
    end
end
t=t-t(start);

measuredTime=t;


%% Extract and plot

[T,X]=extractionReflexHelper(t,q0);
y=X(:,1:4);
for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2))*X(k,3:4)')';
end

figure(100)
clf
subplot(3,1,1)
hold on
plot(xvaf(:,1),xvaf(:,2),'b')
plot(y(:,1),y(:,2),'k')
quiver(y(:,1),y(:,2),xvaf(:,7),xvaf(:,8),'c')
axis equal
subplot(3,1,2)
hold on
vmy=vecmag(y(:,3:4));
plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
subplot(3,1,3)
hold on
plot(t,xvaf(:,2))
quiver(t',xvaf(:,2),xvaf(:,7),xvaf(:,8),'c')

%% Step 1: mark the first peak and somehow optimize straightness up to it
[vals,locs]=findpeaks(vmy,'MINPEAKHEIGHT',.1);

peak=locs(1)-1;

fhigh=find(vecmag(xvaf(:,7:8))>1.5);

subplot(3,1,1)
plot(y(start,1),y(start,2),'go')
plot(y(locs,1),y(locs,2),'rx')
subplot(3,1,2)
plot(t(start),vmy(start),'go')
plot(t(locs),vmy(locs),'rx')
subplot(3,1,3)
plot(t(start),xvaf(start,2),'go')
plot(t(locs),xvaf(locs,2),'rx')
plot(t(fhigh),xvaf(fhigh,2),'r.')

%% First submovement - Use an iterative process to firm up our estimate of the first submovement

%Define a beginning, middle, end, and "direction" range.
%Remember that whole-movement extraction or the remainder from the previous subunit actually lets us go peak-hunting
%Extract with the value of kpgain that we start this iteration with
%Draw up line of launch, 5th order polynomial
%Compare to some notion of convergence. If not converged, adjust kpgain

clear kpg

start=start;
peak_=locs(find(locs>start,1,'first'));
peak=peak_-start;
calib_end=fhigh(1)-start;
range_inds=start:2*peak_-start;
kpg{1}=40*[1 0;0 1];
kpg{2}=.5*[1 0;0 1];
tfull=t;
y1=y;

%close(100)
ITS=2;
while norm(kpg{ITS}-kpg{ITS-1},2)>.01
    kpgain=kpg{ITS};

    t=tfull(range_inds);
    [T,X]=extractionReflexHelper(t,q0);
    y=X(:,1:4);
    for k=1:length(T)
        y(k,1:2)=fkin(X(k,1:2));
        y(k,3:4)=(fJ(X(k,1:2))*X(k,3:4)')';
    end

    figure(ITS)
    set(gcf,'Position', [ 634   430   560   420])
    clf
    subplot(2,1,1)
    hold on
    plot(xvaf(:,1),xvaf(:,2),'b')
    plot(y(:,1),y(:,2),'k')
    quiver(y(:,1),y(:,2),xvaf(range_inds,7),xvaf(range_inds,8),'c')
    plot(y(1,1),y(1,2),'go')
    plot(y(peak,1),y(peak,2),'rx')
    axis equal
    subplot(2,1,2)
    hold on
    vmy=vecmag(y(:,3:4));
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')

    % Make a 5th order polynomial whose span in distance and time is twice the
    % start to the middle. And whose direction is fit from the non-force-perturbed "real" reach.
    fit_inds=1:calib_end;
    line=[y(fit_inds,1) ones(length(fit_inds),1)]\y(fit_inds,2);
    dist=norm((y(1,1:2)-y(peak,1:2))-(dot(y(1,1:2)-y(peak,1:2),[line(1) 1]))*[line(1) 1]);
    tdist=sign(y(peak,1)-y(1,1))*1/(1+line(1)^2)^(1/2)*dist;
    full_inds=1:2*peak;
    [trash,locs_]=findpeaks(vecmag2(y(:,3:4)));
    peak=locs_(1);
    %[trash,peak]=min(abs(t-.15));

    coeff=calcminjerk(y(1,1:2),y(1,1:2)+[2*tdist line(1)*2*tdist],[0 0],[0 0],[0 0],[0 0],0,2*t(peak));
    tmj=t-t(1);
    tmj(tmj>2*t(peak))=2*t(peak);
    [xc,vc,ac]=minjerk(coeff,tmj);
    xc=xc';
    vc=vc';
    ac=ac';

    subplot(2,1,1)
    plot(y(1:peak,1),line(1)*y(1:peak,1)+line(2),'g.')
    plot(y(1,1)+[0 2*tdist],y(1,2)+[0 line(1)*2*tdist],'m')
    plot(xc(:,1),xc(:,2),'m')
    subplot(2,1,2)
    plot(t',vecmag(vc),'m',t,vc(:,1)','m-.',t',vc(:,2),'m--')
    plot(tdist,norm(y(peak,3:4)),'gx')

    tff_launch=xc;
    
    idealVals=zeros(length(t),4);
    
    for k=1:length(tmj)
        q=ikin(xc(k,:));
        qdot=fJq\vc(k,:)';
        qddot=getAlpha(q,qdot,ac(k,:)');
        idealVals(k,:)=[q' qdot'];
        [D,C]=computeDC(q,qdot);
        tff_launch(k,:)=(D*qddot+C)';
    end
    
    errorVals=measuredVals(range_inds,1:4)-idealVals;
    errorTime=measuredTime(range_inds);
    
    tff_ideal=tff_launch;
    tfb_ideal=tff_ideal;
    
    for k=1:length(range_inds)
        [dqi,tfb_ideal(k,:),tff_ideal(k,:)]=armdynamicsInvertedBurdetReflexes(tfull(range_inds(k)),idealVals(k,:)');
    end

    kpgain=kpgain-.1*tfb_ideal(1:peak,:)\(tff_launch(1:peak,:)-tff_ideal(1:peak,:));
    
    
    [T,X]=extractionReflexHelper(t,q0);
    ycor=X(:,1:4);
    for k=1:length(T)
        ycor(k,1:2)=fkin(X(k,1:2));
        ycor(k,3:4)=(fJ(X(k,1:2))*X(k,3:4)')';
    end
    
    subplot(2,1,1)
    plot(ycor(:,1),ycor(:,2),'r')
    subplot(2,1,2)
    plot(t,vecmag(ycor(:,3:4)),'r',t,ycor(:,3),'r-.',t,ycor(:,4),'r--')
    
    %     figure(40+ITS)
    %     clf
    %     hold on
    %     plot(vecmag(tfb),vecmag(tff_launch-tff),'.')
    kpg_scalar=vecmag(tfb(1:peak,:))\vecmag(tff_launch(1:peak,:)-tff(1:peak,:));
    kpg_mat=tfb(1:peak,:)\(tff_launch(1:peak,:)-tff(1:peak,:));
    %     wuh=vecmag(tfb*kpg_mat);
    %     plot(vecmag(tfb),kpg_scalar*vecmag(tfb),'r')
    %     plot(vecmag(tfb),wuh,'r.')
    %
    figure(100)
    set(gcf,'Position',[  1413         339         560         420])
    subplot(3,1,1)
    hold on
    plot(y(:,1),y(:,2),'g')
    axis equal

    kpg{ITS+1}=kpgain-.1*real(kpg_scalar');
    kpg{ITS+1}
    ITS=ITS+1;
end



cleanup
