clc
clear all
close all

load lumps.mat
colors=linspecer(8);
figure(1)
msize=10;

for k=1:45
    for kk=1:8
        if length(lumps(k,kk).dC) > 0
            %lumps(k,kk).backward=lumps(k,kk).L2(1:end-1).*lumps(k,kk).Tn2(1:end-1);
            %lumps(k,kk).forward=lumps(k,kk).L2(2:end).*lumps(k,kk).Tn2(2:end);
            lumps(k,kk).backward=lumps(k,kk).Tn2(1:end-1).^-2;
            lumps(k,kk).forward=lumps(k,kk).Tn2(2:end).^-2;
        end
    end
end

subplot(1,2,1)
hold on
for k=1:8
    subs(k).forward=vertcat(lumps(:,k).forward);
    subs(k).dC=abs([lumps(:,k).dC]').^2;
    plot(subs(k).dC,subs(k).forward,'.','color',colors(k,:),'markersize',msize)
end
title('Forward')
ylabel('Submotion Duration^{-2}')
xlabel('Peak-to-peak Interval squared')

subplot(1,2,2)
hold on
for k=1:8
    subs(k).backward=vertcat(lumps(:,k).backward);
    subs(k).dC=abs([lumps(:,k).dC]').^2;
    plot(subs(k).dC,subs(k).backward,'.','color',colors(k,:),'markersize',msize)
end
title('Backward')
xlabel('Peak-to-peak Interval squared')
