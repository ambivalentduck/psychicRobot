function [y,lumps,L2I1,kpgains]=stiffnessExtract(t,xvaf,N)

global kpgain

%Overestimation of stiffness
y=extract2(t,xvaf,'reflex',[0 inf],[1.5]);
%Speed domain is actually nearly unaffected, but there will be lots of
%residuals
L=findLumps(t,y,length(t)-20,max(findpeaks(vecmag(y)))/3,0);

DIV=3;
while length(L)<2
    DIV=DIV+1;
    L=findLumps(t,y,length(t)-20,max(findpeaks(vecmag(y)))/DIV,0);    
end
L2I1=L(2).inds(1); %To force 2, can just reduce minpeak...should never be a problem

if N>0
    figure(N)
    clf
    hold on
    plot(xvaf(:,1),xvaf(:,2),'k')
    axis equal
    plot(y(:,1),y(:,2),'c')
    for k=1:length(L)
        text(y(L(k).inds(1),1)+k/500,y(L(k).inds(1),2),num2str(k))
        text(y(L(k).inds(end),1)+k/500,y(L(k).inds(end),2),num2str(k))
    end
end

gains=[.2 .3 .5 .7 .9 1.2 1.5];

best=inf;
for G=1:length(gains)
    kpgain=gains(G);
    y=extract(t(1:L2I1),xvaf,'reflex');
    if abs(y(end,2)-.5)<best
        best=abs(y(end,2)-.5);
        besty=y;
        bestG=kpgain;
    end
end

if N>0
    plot(besty(:,1),besty(:,2),'r.')
end

% LMIN=99;
% for G=1:length(gains)
%     kpgain=gains(G);
% 
%     y=extract(t(L2I1:end),xvaf,'reflex',besty(end,:));
%     lumps=findLumps(t,[besty;y(2:end,:)],length(t)-20,max(findpeaks(vecmag(y)))/4);
%     if length(lumps)<LMIN
%         LMIN=length(lumps);
%         y2=[besty;y(2:end,:)];
%         lumps2=lumps;
%         bestKp=kpgain;
%     end
% end
% y=y2;
% lumps=lumps2;

gains2=[.7 .9 1.2 1.5];
eMin=inf;

% %Time?
% for G=1:length(gains2)
%     kpgain=gains(G);
%     y=extract(t(L2I1:end),xvaf,'reflex',besty(end,:));
%     dists=vecmag([y(:,1)-xvaf(end,1) y(:,2)-xvaf(end,2)])<.01;
%     speeds=vecmag(y(:,3:4));
%     f=find((dists<.01)&(speeds<.05),1,'first');
%     if f<eMin
%         eMin=f;
%         y2=[besty;y(2:end,:)];
%         bestKp=kpgain;
%     end
% end

%Pathlength?
for G=1:length(gains2)
    kpgain=gains2(G);
    y=extract(t(L2I1:end),xvaf,'reflex',besty(end,:));
    PL=sum(vecmag(y(2:end,1:2)-y(1:end-1,1:2)));
    dists=vecmag([y(:,1)-xvaf(end,1) y(:,2)-xvaf(end,2)])<.01;
    
    if (PL<eMin)&&~isempty(find(dists))
        eMin=PL;
        y2=[besty;y(2:end,:)];
        bestKp=kpgain;
    end
end

y=y2;
lumps=findLumps(t,y,length(t)-20,max(findpeaks(vecmag(y)))/4);


%findLumps(t,y,length(t)-20,max(findpeaks(vecmag(y)))/4,10);

if N>0
    [C,owned]=lumps2rgbk(lumps);
    for cc=find(owned)'
        plot(y(cc:cc+1,1),y(cc:cc+1,2),'-','Color',C(cc,:))
    end
end

kpgains=[bestG bestKp];