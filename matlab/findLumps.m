function lumps=findLumps(t,y,upperlim,doPlots)

if nargin==2
    upperlim=length(t);
    doPlots=0;
end

if nargin==3
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
ylast=y;

vmy=vecmag(y(:,3:4));
[trash,peak]=max(vmy(1:upperlim));
vals=findpeaks(vmy(1:upperlim),'minpeakheight',.1);

expectedsubs=3*length(vals);
maxlumps=max(10,doPlots);

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

    %Fit a line from beginning to peak, first pass. Dodge issues with
    %vertical lines having inf slope by regressing against time.
    inddist=ceil(.1/gT);
    fit_inds=max(1,peak-inddist):min(length(t),peak+inddist);
    line1=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,1);
    line2=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,2);

    dr=[line1(1) line2(1)];
    dr=dr/norm(dr);
    dp1=zeros(lT,1);
    dpunit=dp1;
    for kk=1:lT
        dp1(kk)=dot(y(kk,3:4),dr);
        dpunit(kk)=dp1(kk)/norm(y(kk,3:4));
    end

    % Second pass...
    gu=gradient(dpunit);
    [vals,locs1]=findpeaks(gu,'minpeakheight',.01);
    [vals,locs2]=findpeaks(-gu,'minpeakheight',.01);

    lower=locs1(find(locs1<peak,1,'last'));
    if isempty(lower)
        lower=1;
    end

    upper=locs2(find(locs2>peak,1,'first'));
    if isempty(upper)
        upper=lT;
    end
 
    dpunit=dpunit; %*norm(y(peak,3:4));
    gu=norm(y(peak,3:4))*gu/max(abs(gu));

    if doPlots
        figure(doPlots+2)
        subplot(ceil((expectedsubs+1)/2),2,k)
        plot(t,dpunit,'.','Color',colors(k,:))
        plot(t,gu,'x','Color',colors(k,:))
        plot(t(lower),norm(y(lower,3:4)),'mo','Markersize',10)
        plot(t(upper),norm(y(upper,3:4)),'m^','Markersize',10)
    end

    unitslope=[-line2(1) line1(1)]/norm([-line2(1) line1(1)]);

    vec1=(y(peak,1:2)-y(lower,1:2))-dot(y(peak,1:2)-y(lower,1:2),unitslope)*unitslope;
    vec2=(y(peak,1:2)-y(upper,1:2))-dot(y(peak,1:2)-y(upper,1:2),unitslope)*unitslope;
    if norm(vec1)<norm(vec2)
        vec=vec1;
        ispan=peak-lower;
    else
        vec=-vec2;
        ispan=upper-peak;
    end
    
    inds=max(1,peak-ispan):min(lT,peak+ispan);
    tspan=t(inds(end))-t(inds(1));
    vmaxreal=norm(y(peak,3:4));
    vmaxcalc=2*norm(vec)*1.875/tspan;
    %if vmaxcalc>vmaxreal
        vec=vec*vmaxreal/vmaxcalc;
    %end

    dp=zeros(length(inds),1);
    for kk=1:length(inds)
        dp(kk)=dot(y(inds(kk),3:4),dr);
    end

    %From line to 5th order polynomial
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
        plot(xc(:,1),xc(:,2),'Color',colors(k,:))
        plot(xc(1,1),xc(1,2),'o','Color',colors(k,:))
        plot(ylast(inds(1),1),ylast(inds(1),2),'x','Color',colors(k,:))
        plot(ylast(inds(end),1),ylast(inds(end),2),'x','Color',colors(k,:))

        subplot(ceil((expectedsubs+1)/2),2,k+1)
        hold on
        plot(y(:,1),y(:,2),'k.')
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
    colormap(lumps2rgbk(lumps));  % This is the color, vary with x in this case.
    col=(1:length(yraw(:,1)))/length(yraw(:,1));
    surface([x';x'],[y';y'],[z';z'],[col;col],'facecol','no','edgecol','interp','linew',2);
    axis equal
end

peaks=[lumps.peak];
[s,i]=sort(peaks);
lumps=lumps(i);
