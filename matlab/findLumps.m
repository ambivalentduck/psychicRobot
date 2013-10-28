function lumps=findLumps(t,y,doPlots)

if nargin<3
    doPlots=0;
end

y(:,1)=y(:,1)-y(1,1);
y(:,2)=y(:,2)-y(1,2);
yraw=y;

vmy=vecmag(y(:,3:4));

%Identify onset of movement
calib=find(vmy<.027); %Square roots are unnecessarily slow.
dtimec=[0 diff(t(calib))]; %Time distance since calibration was satisfied
breaks=find(dtimec>.1);
if isempty(calib)
    start=1;
else
    if isempty(breaks)
        start=calib(end)+1;
    else
        start=calib(breaks(1)-1);
    end
end

[vals,locs]=findpeaks(vmy,'MINPEAKHEIGHT',.1);
locs=locs((y(locs,1)/y(end,1))<.9); 

[trash,start]=min(vmy(1:locs(1)));
t(start)
t(locs(1))

start0=start;
locs0=locs;
y=yraw;
expectedsubs=length(locs)+1;
peak=locs(find(locs>start,1,'first'));

if doPlots
    colors=rand(10*length(locs),3);

    figure(doPlots)
    clf
    subplot(3,1,1)
    hold on
    plot(y(:,1),y(:,2),'k')
    plot(y(start,1),y(start,2),'go')
    plot(y(locs,1),y(locs,2),'rx')
    axis equal
    xlabel('Robot X, m')
    ylabel('Robot Y, m')

    subplot(3,1,2)
    hold on
    vmy=vecmag(y(:,3:4));
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(start),vmy(start),'go')
    plot(t(locs),vmy(locs),'rx')
    xlabel('Time post-onset, sec')
    ylabel('Velocity, m/s')

    figure(doPlots+1)
    clf
    subplot(length(locs)+1,1,1)
    hold on
    plot(y(:,1),y(:,2),'b.')
    plot(y(locs,1),y(locs,2),'rx')

    figure(doPlots+2)
    clf
    subplot(length(locs)+1,1,1)
    hold on
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(start),vmy(start),'go')
    plot(t(locs),vmy(locs),'rx')
end

offset=[0 0];

k=0;
while ((yraw(start,1)/yraw(end,1))<.8)|(k<2)
    k=k+1;
    fit_inds=floor((start+peak)/2):peak;
    full_inds=start:2*peak-start;
    lumps(k).inds=full_inds;
    lumps(k).peak=peak;
    lumps(k).offset=offset;

    %Fit a line from beginning to peak
    line=[y(fit_inds,1) ones(length(fit_inds),1)]\y(fit_inds,2);
    dist=norm((y(start,1:2)-y(peak,1:2))-(dot(y(start,1:2)-y(peak,1:2),[-line(1) 1]))*[-line(1) 1]);
    tdist=sign(y(peak,1)-y(start,1))*1/(1+line(1)^2)^(1/2)*dist;
    
    dr=[tdist line(1)*tdist];
    dr=dr/norm(dr);
    dp=zeros(length(full_inds),1);
    for kk=1:length(full_inds)
        dp(kk)=dot(y(full_inds(kk),3:4),dr);
    end

    %From line to 5th order polynomial
    coeff=calcminjerk(y(start,1:2),y(start,1:2)+[2*tdist line(1)*2*tdist],[0 0],[0 0],[0 0],[0 0],0,t(2*peak-start)-t(start));
    tmj=t(full_inds)-t(fit_inds(1));
    tmj(tmj>2*t(fit_inds(end)))=2*t(fit_inds(end));
    [xc,vc,ac]=minjerk(coeff,tmj);
    xc=xc';
    vc=vc';
    ac=ac';
    lumps(k).y=[xc vc ac];

    if doPlots
        figure(doPlots)
        subplot(3,1,1)
        plot(yraw(fit_inds(1),1)+[0 2*tdist],yraw(fit_inds(1),2)+[0 2*line(1)*tdist],'Color',colors(k,:))
        plot(offset(1)+xc(:,1),offset(2)+xc(:,2),'.','Color',colors(k,:))

        subplot(3,1,2)
        plot(t(full_inds),dp,'.','Color',colors(k,:))
        plot(t(full_inds)',vecmag(vc(:,1:2)),'Color',colors(k,:))
        plot(t(full_inds)',vc(:,1),'-.','Color',colors(k,:))
        plot(t(full_inds)',vc(:,2),'--','Color',colors(k,:))
    end

    %Subtract off 5th order poly
    y(full_inds,:)=y(full_inds,:)-[xc vc];
    y(full_inds(end)+1:end,1:2)=y(full_inds(end)+1:end,1:2)-(xc(end,:)'*ones(1,length(t(full_inds(end)+1:end))))';

    %New start is last time point where x <= 0
    f=find((y(1:2*peak-start,3)<0)&(y(1:2*peak-start,1)<0),1,'last');
    vmy=vecmag(y(:,3:4));
    [vals,locs]=findpeaks(vmy,'MINPEAKHEIGHT',.1);
    locs=locs((y(locs,1)/y(end,1))<.9);

    if doPlots
        figure(doPlots+1)
        subplot(expectedsubs+1,1,k+1)
        hold on
        plot(y(:,1),y(:,2),'.')
        plot(y(locs0,1),y(locs0,2),'mx')
        plot(y(locs,1),y(locs,2),'rx')
        plot(y(f,1),y(f,2),'gx')
        axis equal

        figure(doPlots+2)
        subplot(expectedsubs+1,1,k+1)
        hold on
        plot(t,vmy,'k',t,y(:,3),'k-.',t,y(:,4),'k--')
        plot(t(start),vmy(start),'go')
        plot(t(locs0),vmy(locs0),'mx')
        plot(t(locs),vmy(locs),'rx')
        plot(t(f),vmy(f),'gx')
    end

    start=f;
    peak=locs(find(locs>start,1,'first'));
    offset=offset+xc(end,:);

    lumps(k).resid=y;
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
