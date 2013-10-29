function lumps=findLumps2(t,y,doPlots)

if nargin<3
    doPlots=0;
end

t=t-t(1);
if size(t,2)>size(t,1)
    t=t';
end

y(:,1)=y(:,1)-y(1,1);
y(:,2)=y(:,2)-y(1,2);
yraw=y;
ylast=y;

vmy=vecmag(y(:,3:4));

[trash,peak]=max(vmy);
[trash,f]=findpeaks(-vmy(1:peak));
if isempty(f)
    start=1;
else
    start=f(1);
end

start0=start;

expectedsubs=4+1;

if doPlots
    colors=[1 0 0; 0 1 0; 0 0 1; rand(10,3)];

    figure(doPlots)
    clf
    subplot(3,1,1)
    hold on
    plot(y(:,1),y(:,2),'k')
    plot(y(start,1),y(start,2),'mo')
    plot(y(peak,1),y(peak,2),'mx')
    axis equal
    xlabel('Robot X, m')
    ylabel('Robot Y, m')

    subplot(3,1,2)
    hold on
    vmy=vecmag(y(:,3:4));
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(start),vmy(start),'mo')
    plot(t(peak),vmy(peak),'mx')
    xlabel('Time post-onset, sec')
    ylabel('Velocity, m/s')

    figure(doPlots+1)
    clf
    subplot(expectedsubs+1,1,1)
    hold on
    plot(y(:,1),y(:,2),'k.')

    figure(doPlots+2)
    clf
    subplot(expectedsubs+1,1,1)
    hold on
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(start),vmy(start),'go')
end

offset=0*y(:,1:2);

k=0;
while 1
    k=k+1
    
    %Set
    inddist=floor((peak-start)/4);
    fit_inds=max(1,peak-inddist):min(length(t),peak+inddist);
    full_inds=max(1,start):min(length(t),2*peak-start);
    lumps(k).offset=offset;

    %Fit a line from beginning to peak, first pass. Dodge issues with
    %vertical lines having inf slope by regressing against time.
    line1=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,1);
    line2=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,2);

    dr=[line1(1) line2(1)];
    dr=dr/norm(dr);
    dp1=zeros(length(full_inds),1);
    dpunit=dp1;
    for kk=1:length(full_inds)
        dp1(kk)=dot(y(full_inds(kk),3:4),dr);
        dpunit(kk)=dp1(kk)/norm(y(full_inds(kk),3:4));
    end
    
    % Second pass...
    gu=gradient(dpunit);
    mid=floor(length(full_inds)/2);
    [trash,lower]=min(gu);
    [trash,upper]=max(gu);
    if (lower-mid)>(mid-upper)
        inddist=mid-upper;
    else
        inddist=lower-mid;
    end
    old_inds=full_inds;
    
    if doPlots
        figure(doPlots+2)
        subplot(expectedsubs+1,1,k)
        plot(t(old_inds),dpunit,'.','Color',colors(k,:))
        plot(t(old_inds),gradient(dpunit),'x','Color',colors(k,:))
        plot(t(lower+old_inds(1)-1),norm(ylast(lower+old_inds(1)-1,3:4)),'co','Markersize',5)
        plot(t(upper+old_inds(1)-1),norm(ylast(upper+old_inds(1)-1,3:4)),'cx','Markersize',5)
    end
    
    full_inds=max(1,(peak-inddist)):min(length(t),(peak+inddist));
    start=full_inds(1);

    inddist=floor((peak-start)/4);
    fit_inds=max(1,peak-inddist):min(length(t),peak+inddist);
    line1=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,1);
    line2=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,2);

    unitslope=[-line2(1) line1(1)]/norm([-line2(1) line1(1)]);
    vec1=(y(peak,1:2)-y(start,1:2))-dot(y(peak,1:2)-y(start,1:2),unitslope)*unitslope;
    vec2=(y(peak,1:2)-y(full_inds(end),1:2))-dot(y(peak,1:2)-y(full_inds(end),1:2),unitslope)*unitslope;
    if norm(vec1)<norm(vec2)
        vec=vec1;
    else
        vec=-vec2;
    end

    dr=[line1(1) line2(1)];
    dr=dr/norm(dr);
    dp=zeros(length(full_inds),1);
    for kk=1:length(full_inds)
        dp(kk)=dot(y(full_inds(kk),3:4),dr);
    end
    lumps(k).ownership=zeros(length(t),1);
    lumps(k).ownership(full_inds)=dp;

    %From line to 5th order polynomial
    coeff=calcminjerk(y(peak,1:2)-vec,y(peak,1:2)+vec,[0 0],[0 0],[0 0],[0 0],0,t(full_inds(end))-t(start));
    tmj=t(full_inds)-t(start);
    %tmj(tmj>2*t(fit_inds(end)))=2*t(fit_inds(end));
    [xc,vc,ac]=minjerk(coeff,tmj);
    %Need some logic here to say that if the peak is too high or too low,
    %adjust the time duration since the spacial duration is probably
    %acccurate
    xc=xc';
    vc=vc';
    ac=ac';
    lumps(k).y=[xc vc ac];

    %Subtract off 5th order poly
    ylast=y;
    y(1:full_inds(1)-1,1:2)=[y(1:full_inds(1)-1,1)-xc(1,1) y(1:full_inds(1)-1,2)-xc(1,2)];
    y(full_inds,:)=y(full_inds,:)-[xc vc];
    y(full_inds(end)+1:end,1:2)=y(full_inds(end)+1:end,1:2)-(xc(end,:)'*ones(1,length(t(full_inds(end)+1:end))))';
    vmy=vecmag(y(:,3:4));
    
    if doPlots
        figure(doPlots)
        subplot(3,1,1)
        plot(yraw(peak,1)+[-vec(1) vec(1)],yraw(peak,2)+[-vec(2) vec(2)],'Color',colors(k,:))
        plot(offset(full_inds,1)+xc(:,1),offset(full_inds,2)+xc(:,2),'.','Color',colors(k,:))

        subplot(3,1,2)
        plot(t(old_inds),dp1,'-.','Color',colors(k,:))
        plot(t(full_inds),dp,'.','Color',colors(k,:))
        plot(t(full_inds)',vecmag(vc(:,1:2)),'Color',colors(k,:))
        plot(t(full_inds)',vc(:,1),'-.','Color',colors(k,:))
        plot(t(full_inds)',vc(:,2),'--','Color',colors(k,:))
    
        figure(doPlots+1)
        subplot(expectedsubs+1,1,k)
        plot(ylast(full_inds,1),ylast(full_inds,2),'.','Color',colors(k,:))
        plot(xc(:,1),xc(:,2),'Color',colors(k,:))
        plot(xc(1,1),xc(1,2),'o','Color',colors(k,:))
        plot(ylast(full_inds(1),1),ylast(full_inds(1),2),'x','Color',colors(k,:))
        plot(ylast(full_inds(end),1),ylast(full_inds(end),2),'x','Color',colors(k,:))

        subplot(expectedsubs+1,1,k+1)
        hold on
        plot(y(:,1),y(:,2),'k.')
        plot(y(1,1),y(1,2),'mo')
        axis equal

        figure(doPlots+2)
        subplot(expectedsubs+1,1,k)
        plot(t(full_inds),vecmag(vc),'Color',colors(k,:))

        subplot(expectedsubs+1,1,k+1)
        hold on
        plot(t,vmy,'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    end

    lumps(k).inds=full_inds;
    lumps(k).peak=peak;

    [pkheight,peak]=max(vmy);
    [trash,f]=findpeaks(-vmy(1:peak));
    if isempty(f)
        start=1;
    else
        start=f(end);
    end
    if doPlots
        figure(doPlots+2)
    
        subplot(expectedsubs+1,1,k+1)
        hold on
        plot(t(start),vmy(start),'mo')
        plot(t(peak),vmy(peak),'mx')
    end
    

    offset(1:full_inds(1)-1,:)=[offset(1:full_inds(1)-1,1)+xc(1,1) offset(1:full_inds(1)-1,2)+xc(1,2)];
    offset(full_inds,:)=offset(full_inds,:)+xc;
    offset(full_inds(end)+1:end,:)=[offset(full_inds(end)+1:end,1)+xc(end,1) offset(full_inds(end)+1:end,2)+xc(end,2)];

    lumps(k).resid=y;
    
    if pkheight<.1
        break
    end
end

if doPlots
    figure(doPlots)
    subplot(3,1,3)
    x = yraw(:,1);
    y = yraw(:,2);
    z = zeros(size(x));
    colormap(lumps2rgbk(lumps,yraw));  % This is the color, vary with x in this case.
    col=(1:length(yraw(:,1)))/length(yraw(:,1));
    surface([x';x'],[y';y'],[z';z'],[col;col],'facecol','no','edgecol','interp','linew',2);
end
