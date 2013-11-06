clc
clear all

load('exampledecomp');

upperlim=mark-10;
t=t(1:115);
y=y(1:115,:);

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

colors=colorScheme(10);

peakOrder=[2 4 1 3];

figure(1)
clf

for k=1:6
    subplot(6,2,2*k)
    plot([0 1],[0 0],'w','linewidth',2)
end

subplot(6,2,1)
hold on
%plot(y(:,1),y(:,2),'k')
axis equal

subplot(6,2,2)
hold on
%vmy=vecmag(y(:,3:4));
%plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')

subplot(6,2,3)
hold on
plot(y(:,1),y(:,2),'k.')
plot(y(peak,1),y(peak,2),'mx','Markersize',10)
axis equal

subplot(6,2,4)
hold on
plot(t,vecmag(y(:,3:4)),'k',t,y(:,3),'k-.',t,y(:,4),'k--')
plot(t(peak),vecmag(y(peak,3:4)),'mx','Markersize',10)

for k=1:4
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
    [vals,locs1]=findpeaks(gu1,'minpeakheight',.01);
    [vals,locs2]=findpeaks(-gu1,'minpeakheight',.01);
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

    unitslope=[-line2(1) line1(1)]/norm([-line2(1) line1(1)]);
    vec1=(y(peak,1:2)-y(lower2,1:2))-dot(y(peak,1:2)-y(lower2,1:2),unitslope)*unitslope;
    vec2=(y(peak,1:2)-y(upper2,1:2))-dot(y(peak,1:2)-y(upper2,1:2),unitslope)*unitslope;
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

    tspan2=tspan*vmaxcalc/vmaxreal;
    ispan2=round(tspan2/(2*gT));
    lower3=max(1,peak-ispan2);
    upper3=min(lT,peak+ispan2);
    inds=lower3:upper3;
    tspan=t(inds(end))-t(inds(1));

    dpunit1=dpunit1*norm(y(peak,3:4));

    subplot(6,2,(k+1)*2)
    plot(t,dpunit1,'--','Color',colors(peakOrder(k),:))
    plot(t,gu2,'.','Color',colors(peakOrder(k),:))
    plot(t(lower2),norm(y(lower2,3:4)),'m^','Markersize',10)
    plot(t(upper2),norm(y(upper2,3:4)),'m^','Markersize',10)
    plot(t(lower3),norm(y(lower3,3:4)),'mo','Markersize',10)
    plot(t(upper3),norm(y(upper3,3:4)),'mo','Markersize',10)

    %Third pass
    xi=(y(peak,1:2)-vec)';
    xf=(y(peak,1:2)+vec)';
    ta=(t(inds)'-t(inds(1)))/tspan;
    xc=(xi*ones(size(ta))+(xf-xi)*(10*ta.^3-15*ta.^4+6*ta.^5))';
    vc=((xf-xi)*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';
    ac=((xf-xi)*(60*ta-180*ta.^2+120*ta.^3)/(tspan^2))';
    lumps(k).y=[xc vc ac];

    lumps(k).ownership=zeros(length(t),1);
    lumps(k).ownership(inds)=(30*ta.^2-60*ta.^3+30*ta.^4)/1.875;

    %Subtract off 5th order poly
    ylast=y;
    y(1:inds(1)-1,1:2)=[y(1:inds(1)-1,1)-xc(1,1) y(1:inds(1)-1,2)-xc(1,2)];
    y(inds,:)=y(inds,:)-[xc vc];
    y(inds(end)+1:end,1:2)=[y(inds(end)+1:end,1)-xc(end,1) y(inds(end)+1:end,2)-xc(end,2)];
    vmy=vecmag(y(:,3:4));

    subplot(6,2,k*2+1)
    plot(ylast(peak,1)+[-vec(1) vec(1)],ylast(peak,2)+[-vec(2) vec(2)],'Color',colors(peakOrder(k),:),'LineWidth',3)
    
    subplot(6,2,(k+1)*2)
    plot(t(inds)',vecmag(vc(:,1:2)),'Color',colors(peakOrder(k),:))
    plot(t(inds)',vc(:,1),'-.','Color',colors(peakOrder(k),:))
    plot(t(inds)',vc(:,2),'--','Color',colors(peakOrder(k),:))

    subplot(6,2,k*2+3)
    hold on
    plot(y(:,1),y(:,2),'k.')
    axis equal

    subplot(6,2,k*2+4)
    hold on
    plot(t,vmy,'k',t,y(:,3),'k-.',t,y(:,4),'k--')

    lumps(k).inds=inds;
    lumps(k).peak=peak;

    [pkheight,peak]=max(vmy(1:upperlim));

    subplot(6,2,k*2+3)
    plot(y(peak,1),y(peak,2),'mx','Markersize',10)
    subplot(6,2,k*2+4)
    plot(t(peak),vmy(peak),'mx','Markersize',10)


    lumps(k).resid=y;

    if pkheight<.1
        break
    end
end

peaks=[lumps.peak];
[s,i]=sort(peaks);
lumps=lumps(i);

subplot(6,2,1)
hold on
subplot(6,2,2)
hold on


y = yraw;
vmy=vecmag(yraw(:,3:4));
[C,owned]=lumps2rgbk(lumps);
owned=find(owned(owned~=0)');
for cc=owned
    subplot(6,2,1)
    plot(y(cc:cc+1,1),y(cc:cc+1,2),'-','Color',C(cc,:),'linewidth',2)
    subplot(6,2,2)
    plot(t(cc:cc+1),vmy(cc:cc+1),'-','Color',C(cc,:),'linewidth',2)
end

subplot(6,2,1)
axis equal
for k=1:4
    subplot(6,2,1)
    offset=yraw(lumps(k).peak,1:2)-lumps(k).y(floor(end/2)+1,1:2);
    plot(lumps(k).y(:,1)+offset(1),lumps(k).y(:,2)+offset(2),'color',colors(k,:),'linewidth',3)
    subplot(6,2,2)
    plot(t(lumps(k).inds),vecmag(lumps(k).y(:,3:4)),'color',colors(k,:),'linewidth',3)
end


xl=[-.05 .3];
yl=[-.1 .01];
for k=1:6
    subplot(6,2,-1+k*2)
    xlim(xl);
    ylim(yl);
    if k==1
        plot([0 0],[0 -.05],'w','linewidth',3)
        text(-.01,-.04,'5 cm','color','w','rotation',90)
    end
    axis off
    subplot(6,2,k*2)
    xlim([-.1 1.2]);
    ylim([-.5 .5]);
    axis off
    if k==1
        text(.95,-.09,'1s','color','w')
        plot([0 0],[0 .5],'w','linewidth',2)
        text(-.03,0,'50 cm/s','color','w','rotation',90)
        
    end
end
backgray=.6;
set(gcf,'Color',backgray*[1 1 1])
