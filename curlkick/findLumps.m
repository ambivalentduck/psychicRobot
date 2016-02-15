function [lumps,y]=findLumps(t,y,searchinds)

%Accept inputs of the form:
% column t
% columns y rotated to be progress (0 to Y), error
% searchinds lets you pad and not worry about finding lumps outside the reach

doPlots=0;
gT=mean(gradient(t));
lT=length(t);
trueinds=1:lT;

%Make plotting easier for later
y(:,1)=y(:,1)-y(1,1);
y(:,2)=y(:,2)-y(1,2);
yraw=y;

colors=[1 0 0;
    0 1 0;
    0 0 1;
    .8 .8 0;
    0 1 1;
    1 0 1;
    rand(10,3)];

if doPlots
    figure(1)
    clf
    hold on
    vmy=vecmag(y);
    plot(t(searchinds),vmy(searchinds),'k',t,y(:,1),'k-.',t,y(:,2),'k--')
    xlabel('Time post-onset, sec')
    ylabel('Velocity, m/s')
end

k=0;
minpeakheight=.03;
mingradpeakheight=.01;

while k<5
    k=k+1;
    
    %Assume the biggest peak is always pure
    % v_peak=L*1.875/S
    % L is rough, but S we can get from a dot product.
    vmy=vecmag(y(searchinds,:));
    [val,peak]=max(vmy);
    peak=trueinds(searchinds(peak));
    %Make sure we're not catching a rising or falling peak from the
    %previous movement
    %     if ~any(peak==searchinds([1 end]))
    %         if val<minpeakheight
    %             break
    %         end
    %     else
    %         [vals,peaks]=findpeaks(vmy,'minpeakheight',minpeakheight);
    %         if isempty(peaks)
    %             break
    %         end
    %         kk=1;
    %         [svals,si]=sort(vals);
    %         peaks=peaks(si);
    %         while kk<length(vals)
    %
    %
    %
    %     end
    
    %Shitty coding, grossly inefficient
    dr=y(peak,:);
    dr=dr/norm(dr);
    dpunit=zeros(lT,1);
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
    %S should be the lower of 2(upper-C) or 2(C-lower)
    C=t(peak);
    %C=(t(upper)+t(lower))/2;
    %S=gT*(upper-lower)*1.05; %Probably add a fudge factor
    S=min(max(2*min(t(upper)-C,C-t(lower))*1.1,gT),1);
    L=y(peak,:)*S/1.875;
    lumps(k).C=C;
    lumps(k).S=S;
    lumps(k).L=L;
    
    tau=(t-C)/S+.5;
    tau=max(min(tau,1),0);
    kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/S;
    y=y-kappa*L;
    
    if doPlots
        plot(t(searchinds),kappa(searchinds)*L,'-','color',colors(k,:))
        plot(t(searchinds),vecmag(y(searchinds,:)),'-.','color',colors(k,:))
    end
end

if doPlots
    plot(t,vecmag(y),'c-')
end
