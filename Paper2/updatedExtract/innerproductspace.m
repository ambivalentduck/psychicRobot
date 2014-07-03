clc
clear all
close all

%% Set them up

%Presume the shortest human submovement is 150ms and the longest is 450ms.
tspan=[.15 .45];

%Presume minima and maxima in speed, too
mspan=[.01 .05];

Nsubs=20;
ts=zeros(Nsubs,2);
dirs=zeros(Nsubs,2);
mags=zeros(Nsubs,1);

for k=1:Nsubs
    if k==1
        ts(k,1)=0;
    elseif k==2
        ts(k,1)=rand*(tspan(2)-tspan(1))+tspan(1);
    else
        if rand>.3
            ts(k,1)=ts(k-2,2);
        else
            ts(k,1)=ts(k-2,2)+rand*tspan(1);
        end
    end
    ts(k,2)=ts(k,1)+rand*(tspan(2)-tspan(1))+tspan(1);
    
    mags(k)=rand*(mspan(2)-mspan(1))+mspan(1);
    dirs(k,:)=2*rand(1,2)-1;
    dirs(k,:)=dirs(k,:)/norm(dirs(k,:));
end

figure(1)
clf
subplot(3,1,1)
hold on
plot(ts(:,1),zeros(Nsubs,1),'bo',ts(:,2),zeros(Nsubs,1),'bx')

%% Get a speed profile.

tstep=.001;
t=(0:tstep:max(max(ts)))';
v=zeros(length(t),2);

for k=1:Nsubs
    fLow=find(t>=ts(k,1),1,'first');
    if isempty(fLow)
        fLow=1;
    end
    fHigh=find(t<=ts(k,2),1,'last');
    if isempty(fHigh)
        fHigh=length(t);
    end
    
    ta=(t(fLow:fHigh)-ts(k,1))/(ts(k,2)-ts(k,1));
    vd=(30*ta.^2-60*ta.^3+30*ta.^4)*(mags(k)/(ts(k,2)-ts(k,1)))*dirs(k,:);
    c=.9*rand(3,1)+.1;
    plot(t(fLow:fHigh),sqrt(sum(vd.^2,2)),'-','color',c,'linewidth',2)
    v(fLow:fHigh,:)=v(fLow:fHigh,:)+vd;
end

plot(t,sqrt(sum(v.^2,2)),'k')

%% Get a position profile

x=cumsum([0 .5;v*tstep],1);
figure(2)
clf
plot(x(:,1),x(:,2))

%% Attempt FIFO decomposition of the velocity profile.

%while there's data left to fit
c=0;
lookahead=round(.05/tstep);
dotspan=round(.25/tstep);

figure(1)
subplot(3,1,3)
hold on
L=1;

subplot(3,1,2)
hold on

vn=v;

while c<length(t)
    c=c+1;
    %fit 50 ms
    m=vn(c+lookahead-1,:);
    m=m/norm(m);
    D=min(length(t),c+dotspan);
    dotspan=sum(vn(c:D,1)*m(1)+vn(c:D,2)*m(2),2)/norm(m);
    subplot(3,1,3)
    rcolor=rand(1,3);
    plot(t(c:D),dotspan,'color',rcolor)
    [val,loc]=max(dotspan);
    plot(t(c+loc-1),val,'rx')
    xlim([0 t(end)])
    drawnow
    subplot(3,1,2)
    lump(L).mag=val;
    lump(L).dir=m;
    lump(L).ts=[t(c) t(c-1+2*loc)];
    vn(c:c-1+loc,:)=vn(c:c-1+loc,:)-dotspan(1:loc)*m;
    vn(c+loc:c-1+2*loc,:)=vn(c+loc:c-1+2*loc,:)-wrev(dotspan(1:loc))*m;
    plot(t,sqrt(sum(vn.^2,2)),'color',rcolor)
    L=L+1
    c=find(sum(vn.^2,2)>.01,1,'first');
end
