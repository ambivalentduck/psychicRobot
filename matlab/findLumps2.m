function lumps=findLumps2(t,y,upperlim,doPlots)

if nargin==2
    upperlim=length(t);
    doPlots=0;
end

if nargin==3
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

[trash,peak]=max(vmy(1:upperlim));
[trash,f]=findpeaks(-vmy(1:peak));
if isempty(f)
    start=1;
else
    start=f(1);
end
inddist=peak-start+1;

expectedsubs=min(10,doPlots);
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
    axis equal

    figure(doPlots+2)
    clf
    subplot(expectedsubs+1,1,1)
    hold on
    plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
    plot(t(start),vmy(start),'go')
end

offset=0*y(:,1:2);

k=0;
while k<maxlumps
    k=k+1;

    %Set
    full_inds=max(1,(peak-inddist)):min(length(t),(peak+inddist));
    inddist=floor((peak-start)/4);
    fit_inds=max(1,peak-inddist):min(length(t),peak+inddist);
    
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
    [lval,lower]=min(gu(mid:end));
    lower=lower+mid-1;
    [uval,upper]=max(gu(1:mid));
    
    gradthresh=.01;
    if (abs(lval)<gradthresh)&&(abs(uval)<gradthresh)
        figure_out_dynamic_expansion=1
        lval
        uval
        return
    end
    
    if abs(lval)<gradthresh
        inddist=mid-upper;
    elseif abs(uval)<gradthresh
        inddist=lower-mid;
    elseif (lower-mid)>(mid-upper)
        inddist=mid-upper;
    else
        inddist=lower-mid;
    end
    
    old_inds=full_inds;
    dpunit=dpunit*norm(y(peak,3:4));
    gu=norm(y(peak,3:4))*gu/max(abs(gu));

    if doPlots
        figure(doPlots+2)
        subplot(expectedsubs+1,1,k)
        plot(t(old_inds),dpunit,'.','Color',colors(k,:))
        plot(t(old_inds),gu,'x','Color',colors(k,:))
        plot(t(lower+old_inds(1)-1),norm(y(lower+old_inds(1)-1,3:4)),'co','Markersize',10)
        plot(t(upper+old_inds(1)-1),norm(y(upper+old_inds(1)-1,3:4)),'cx','Markersize',10)
    end

    start=max(1,(peak-inddist));

    fitinddist=floor((peak-start)/4);
    fit_inds=max(1,peak-fitinddist):min(length(t),peak+fitinddist);
    line1=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,1);
    line2=[t(fit_inds) ones(length(fit_inds),1)]\y(fit_inds,2);
    unitslope=[-line2(1) line1(1)]/norm([-line2(1) line1(1)]);

    full_inds=max(1,(peak-inddist)):min(length(t),(peak+inddist));

    vec1=(y(peak,1:2)-y(full_inds(1),1:2))-dot(y(peak,1:2)-y(full_inds(1),1:2),unitslope)*unitslope;
    vec2=(y(peak,1:2)-y(full_inds(end),1:2))-dot(y(peak,1:2)-y(full_inds(end),1:2),unitslope)*unitslope;
    if norm(vec1)<norm(vec2)
        vec=vec1;
    else
        vec=-vec2;
    end
    
    tspan=t(full_inds(end))-t(full_inds(1));
    vmaxreal=norm(y(peak,3:4));
    vmaxcalc=2*norm(vec)*1.875/tspan;
    %if vmaxcalc>vmaxreal
        vec=vec*vmaxreal/vmaxcalc;
    %end
    
    dr=[line1(1) line2(1)];
    dr=dr/norm(dr);
    dp=zeros(length(full_inds),1);
    for kk=1:length(full_inds)
        dp(kk)=dot(y(full_inds(kk),3:4),dr);
    end

    %From line to 5th order polynomial
    xi=(y(peak,1:2)-vec)';
    xf=(y(peak,1:2)+vec)';
    ta=(t(full_inds)'-t(full_inds(1)))/tspan;
    xc=(xi*ones(size(ta))+(xf-xi)*(10*ta.^3-15*ta.^4+6*ta.^5))';
    vc=((xf-xi)*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';
    ac=((xf-xi)*(60*ta-180*ta.^2+120*ta.^3)/(tspan^2))';
    lumps(k).y=[xc vc ac];
    
    lumps(k).ownership=zeros(length(t),1);
    lumps(k).ownership(full_inds)=(30*ta.^2-60*ta.^3+30*ta.^4)/1.875;

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
        %plot(t(full_inds),dp,'.','Color',colors(k,:))
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

    [pkheight,peak]=max(vmy(1:upperlim));
    [trash,f]=findpeaks(-vmy);
    if isempty(f)
        inddist=floor(1.2*min(peak,length(t)-peak));
    else
        lower=f(find(f<peak,1,'last')); %Start is last low before peak
        upper=f(find(f>peak,1,'first')); %Start is last low before peak
        if isempty(lower)
            lower=1;
        end
        if isempty(upper)
            upper=length(t);
        end
        inddist=max(peak-lower,upper-peak); %Back up a tiny bit to make sure you get it.
    end
    start=max(1,peak-inddist+1);
    
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
    colormap(lumps2rgbk(lumps));  % This is the color, vary with x in this case.
    col=(1:length(yraw(:,1)))/length(yraw(:,1));
    surface([x';x'],[y';y'],[z';z'],[col;col],'facecol','no','edgecol','interp','linew',2);
    axis equal
end

peaks=[lumps.peak];
[s,i]=sort(peaks);
lumps=lumps(i);
