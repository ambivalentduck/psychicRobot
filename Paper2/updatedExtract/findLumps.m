function lumps=findLumps(t,y,upperlim,minpeakheight,maxlumps,doPlots)

switch nargin
    case 2
        upperlim=length(t);
        minpeakheight=.1;
        maxlumps=10;
        doPlots=0;
    case 3
        minpeakheight=.1;
        maxlumps=10;
        doPlots=0;
    case 4
        maxlumps=10;
        doPlots=0;
    case 5
        doPlots=0;
end

t=t-t(1);
gT=mean(gradient(t));
lT=length(t);
if size(t,2)>size(t,1)
    t=t';
end

%Make plotting easier for later
y(:,1)=y(:,1)-y(1,1);
y(:,2)=y(:,2)-y(1,2);
yraw=y;

vmy=vecmag(y(:,3:4));
[trash,peak]=max(vmy(1:upperlim));
vals=findpeaks(vmy(1:upperlim),'minpeakheight',minpeakheight);

expectedsubs=maxlumps;

if doPlots
    colors=[1 0 0;
        0 1 0;
        0 0 1;
        .8 .8 0;
        0 1 1;
        1 0 1;
        rand(10,3)];

    figure(doPlots)
    clf
    subplot(3,1,1)
    hold on
    plot(y(:,1),y(:,2),'k')
    plot(y(peak,1),y(peak,2),'mx')
    axis equal
    xlabel('Robot X, m')
    ylabel('Robot Y, m')

    subplot(3,1,2)
    hold on
    vmy=vecmag(y(:,3:4));
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(peak),vmy(peak),'mx')
    xlabel('Time post-onset, sec')
    ylabel('Velocity, m/s')

    figure(doPlots+1)
    clf
    subplot(ceil((expectedsubs+1)/2),2,1)
    hold on
    plot(y(:,1),y(:,2),'k.')
    axis equal

    figure(doPlots+2)
    clf
    subplot(ceil((expectedsubs+1)/2),2,1)
    hold on
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(peak),vecmag(y(peak,3:4)),'mx','Markersize',10)
end

offset=0*y(:,1:2);

k=0;
while k<maxlumps
    k=k+1;

    %Three passes: gradient from peak direction gives line. Line gives min-jerk fit
    %There's a quiet assumption that the found peak is perfect and that
    %getting things symmetrical about it is key. Maybe fix in a v2.0.

    %First pass: get a fitting interval via dot() with peak.
    dr=y(peak,3:4);
    dr=dr/norm(dr);
    dp1=zeros(lT,1);
    dpunit1=dp1;
    for kk=1:lT
        dp1(kk)=dot(y(kk,3:4),dr);
        dpunit1(kk)=dp1(kk)/norm(y(kk,3:4));
    end
    gu1=gradient(dpunit1);
    [vals,locs1]=findpeaks(gu1,'minpeakheight',minpeakheight);
    [vals,locs2]=findpeaks(-gu1,'minpeakheight',minpeakheight);
    lower1=locs1(find(locs1<peak,1,'last'));
    upper1=locs2(find(locs2>peak,1,'first'));
    if isempty(lower1)
        lower1=1;
    end
    if isempty(upper1)
        upper1=lT;
    end
    ispan=ceil(min(peak-lower1,upper1-peak)/2); % /2 is not obvious...optimize.
    fit_inds=max(1,peak-ispan):min(lT,peak+ispan);

    %Second pass: Fit a line to first pass interval. Dodge issues with
    %vertical lines having inf slope by regressing against time.
    line1=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,1);
    line2=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,2);
    dr=[line1(1) line2(1)];
    dr=dr/norm(dr);
    dp2=zeros(lT,1);
    dpunit2=dp2;
    for kk=1:lT
        dp2(kk)=dot(y(kk,3:4),dr);
        dpunit2(kk)=dp2(kk)/norm(y(kk,3:4));
    end
    gu2=gradient(dpunit2);
    [vals,locs1]=findpeaks(gu2,'minpeakheight',.01);
    [vals,locs2]=findpeaks(-gu2,'minpeakheight',.01);
    lower2=locs1(find(locs1<peak,1,'last'));
    upper2=locs2(find(locs2>peak,1,'first'));
    if isempty(lower2)
        lower2=1;
    end
    if isempty(upper2)
        upper2=lT;
    end
    gu2=norm(y(peak,3:4))*gu2/max(abs(gu2));

    vec1=dot(y(peak,1:2)-y(lower2,1:2),dr)*dr;
    vec2=dot(y(peak,1:2)-y(upper2,1:2),dr)*dr;
    if norm(vec1)<norm(vec2)
        vec=vec1;
        ispan=peak-lower2;
    else
        vec=-vec2;
        ispan=upper2-peak;
    end

    inds=max(1,peak-ispan):min(lT,peak+ispan);
    tspan=t(inds(end))-t(inds(1));
    vmaxreal=norm(y(peak,3:4));
    vmaxcalc=2*norm(vec)*1.875/tspan;
    if vmaxcalc>vmaxreal
        tspan2=tspan*vmaxcalc/vmaxreal;
        ispan2=round(tspan2/(2*gT));
        lower3=max(1,peak-ispan2);
        upper3=min(lT,peak+ispan2);
        inds=lower3:upper3;
        tspan=t(inds(end))-t(inds(1));
    else
        lower3=lower2;
        upper3=upper2;
    end

    if doPlots
        figure(doPlots+2)
        subplot(ceil((expectedsubs+1)/2),2,k)
        %plot(t,dpunit1,'--','Color',colors(k,:))
        %plot(t,dpunit2,'-.','Color',colors(k,:))
        %plot(t,gu1,'v','Color',colors(k,:))
        %plot(t,gu2,'^','Color',colors(k,:))
        plot(t(lower1),norm(y(lower1,3:4)),'m<','Markersize',10)
        plot(t(upper1),norm(y(upper1,3:4)),'m>','Markersize',10)
        plot(t(lower2),norm(y(lower2,3:4)),'mo','Markersize',10)
        plot(t(upper2),norm(y(upper2,3:4)),'mx','Markersize',10)
        plot(t(lower3),norm(y(lower3,3:4)),'ko','Markersize',10)
        plot(t(upper3),norm(y(upper3,3:4)),'kx','Markersize',10)
    end

    %Third pass
    xi=(y(peak,1:2)-vec)';
    xf=(y(peak,1:2)+vec)';
    ta=(t(inds)'-t(inds(1)))/tspan;
    xc=(xi*ones(size(ta))+(xf-xi)*(10*ta.^3-15*ta.^4+6*ta.^5))';
    vc=((xf-xi)*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';
    ac=((xf-xi)*(60*ta-180*ta.^2+120*ta.^3)/(tspan^2))';
    lumps(k).y=[xc vc ac];

    lumps(k).offset=offset;
    lumps(k).ownership=zeros(length(t),1);
    lumps(k).ownership(inds)=(30*ta.^2-60*ta.^3+30*ta.^4)/1.875;

    %Subtract off 5th order poly
    ylast=y;
    y(1:inds(1)-1,1:2)=[y(1:inds(1)-1,1)-xc(1,1) y(1:inds(1)-1,2)-xc(1,2)];
    y(inds,:)=y(inds,:)-[xc vc];
    y(inds(end)+1:end,1:2)=[y(inds(end)+1:end,1)-xc(end,1) y(inds(end)+1:end,2)-xc(end,2)];
    vmy=vecmag(y(:,3:4));

    if doPlots
        figure(doPlots)
        subplot(3,1,1)
        plot(yraw(peak,1)+[-vec(1) vec(1)],yraw(peak,2)+[-vec(2) vec(2)],'Color',colors(k,:))
        plot(offset(inds,1)+xc(:,1),offset(inds,2)+xc(:,2),'.','Color',colors(k,:))

        subplot(3,1,2)
        plot(t,dp1,'-.','Color',colors(k,:))
        %plot(t(inds),dp,'.','Color',colors(k,:))
        plot(t(inds)',vecmag(vc(:,1:2)),'Color',colors(k,:))
        plot(t(inds)',vc(:,1),'-.','Color',colors(k,:))
        plot(t(inds)',vc(:,2),'--','Color',colors(k,:))

        figure(doPlots+1)
        subplot(ceil((expectedsubs+1)/2),2,k)
        plot(ylast(inds,1),ylast(inds,2),'.','Color',colors(k,:))
        %plot(xc(:,1),xc(:,2),'Color',colors(k,:),'Linewidth',3)
        plot(xc(ceil(end/2),1),xc(ceil(end/2),2),'x','Color',colors(k,:),'Markersize',10)
        plot(xc(1,1),xc(1,2),'o','Color',colors(k,:))
        plot(ylast(inds(1),1),ylast(inds(1),2),'x','Color',colors(k,:))
        plot(ylast(inds(end),1),ylast(inds(end),2),'x','Color',colors(k,:))

        subplot(ceil((expectedsubs+1)/2),2,k+1)
        hold on
        plot(y(1:upperlim,1),y(1:upperlim,2),'k.')
        plot(y(1,1),y(1,2),'mo')
        axis equal

        figure(doPlots+2)
        subplot(ceil((expectedsubs+1)/2),2,k)
        plot(t(inds),vecmag(vc),'Color',colors(k,:))

        subplot(ceil((expectedsubs+1)/2),2,k+1)
        hold on
        plot(t,vmy,'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    end

    lumps(k).inds=inds;
    lumps(k).peak=peak;

    [pkheight,peak]=max(vmy(1:upperlim));

    if doPlots
        figure(doPlots+2)
        subplot(ceil((expectedsubs+1)/2),2,k+1)
        plot(t(peak),vmy(peak),'mx','Markersize',10)
    end

    offset(1:inds(1)-1,:)=[offset(1:inds(1)-1,1)+xc(1,1) offset(1:inds(1)-1,2)+xc(1,2)];
    offset(inds,:)=offset(inds,:)+xc;
    offset(inds(end)+1:end,:)=[offset(inds(end)+1:end,1)+xc(end,1) offset(inds(end)+1:end,2)+xc(end,2)];

    lumps(k).resid=y;

    if pkheight<minpeakheight
        break
    end
end

if doPlots
    figure(doPlots)
    subplot(3,1,3)
    x = yraw(:,1);
    y = yraw(:,2);
    z = zeros(size(x));
    colormap(lumps2rgbk(lumps));  % This is the color, vary with x in this case.
    col=(1:length(yraw(:,1)))/length(yraw(:,1));
    surface([x';x'],[y';y'],[z';z'],[col;col],'facecol','no','edgecol','interp','linew',2);
    axis equal
end

peaks=[lumps.peak];
[s,i]=sort(peaks);
lumps=lumps(i);
