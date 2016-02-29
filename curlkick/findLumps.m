function [lumps,y]=findLumps(t,y,searchinds,varargin)

%Accept inputs of the form:
% column t
% columns y rotated to be progress (0 to Y), error
% searchinds lets you pad and not worry about finding lumps outside the reach

if nargin<4
    maxlumps=15;
else
    maxlumps=varargin{1};
end

doPlots=0;
gT=mean(gradient(t));
lT=length(t);
trueinds=1:lT;

if doPlots
    colors=[1 0 0;
    0 1 0;
    0 0 1;
    .8 .8 0;
    0 1 1;
    1 0 1;
    rand(10,3)];

    figure(117)
    clf
    hold on
    vmy=vecmag(y);
    plot(t,vmy,'color',[.5 .5 .5])
    plot(t(searchinds),vmy(searchinds),'k',t,y(:,1),'k-.',t,y(:,2),'k--')
    xlabel('Time post-onset, sec')
    ylabel('Velocity, m/s')
end

k=0;
minsublength=.005; %5 mm
mingradpeakheight=.002;

while k<maxlumps
    k=k+1;
    
    %Assume the biggest peak is always pure
    % v_peak=L*1.875/S
    % L is rough, but S we can get from a dot product.
    vmy=vecmag(y(searchinds,:));
    [val,peak]=max(vmy);
    peak=trueinds(searchinds(peak));
    
    
    %Shitty coding, grossly inefficient
    dr=y(peak,:);
    dr=dr/norm(dr);
    for kk=1:lT
        dp(kk)=dot(y(kk,:),dr)/norm(y(kk,:));
    end
    gu=gradient(dp);
    [~,locs1]=findpeaks(gu,'minpeakheight',mingradpeakheight);
    [~,locs2]=findpeaks(-gu,'minpeakheight',mingradpeakheight);
    lower=locs1(find(locs1<peak,1,'last'));
    upper=locs2(find(locs2>peak,1,'first'));
    if isempty(lower)
        lower=1;
    end
    if isempty(upper)
        upper=lT;
    end
    
    C=t(peak);
    S=min(1.3*max(2*min(t(upper)-C,C-t(lower)),gT),1.2);
    %S=1.2*max(2*min(t(upper)-C,C-t(lower)),gT);
    L=y(peak,:)*S/1.875; %1.875 corresponds to min jerk *only*
    
    if (norm(L)<minsublength)&(nargin<4)
        break
    end
    
    lumps(k).C=C;
    lumps(k).S=S;
    lumps(k).L=L;
    
    ylump=lumpy(t,lumps(k));
    y=y-ylump;
    
    if doPlots
        plot(t,vecmag(ylump),'-','color',colors(k,:))
        plot(t,vecmag(y),'-.','color',colors(k,:))
    end
end

if doPlots
    plot(t,vecmag(y),'c-')
end
