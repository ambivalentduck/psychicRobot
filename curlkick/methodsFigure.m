clc
clear all

load ../Data/curlkick/curlkick1Y.mat

f=find(([trials.targetcat]~=0)&([trials.disturbcat]));

inds=[22, 20, 38, 17, 11];
T=f(inds(5));

figure(1)
clf

y=trials(T).y;
t=trials(T).ty;
inds=trials(T).i0:trials(T).if+20;
[~,peaks]=findpeaks(-vecmag(y(inds,3:4)));
inds=inds(1):inds(peaks(end)-5);
t=t-t(inds(1));

[lumps,resid]=findLumps(t',y(:,3:4),inds);

[~,peakinds]=findpeaks(vecmag(y(:,3:4)));

subplot(2,1,1)
hold on
plot(y(:,1),y(:,2),'k')
plot(y(peakinds,1),y(peakinds,2),'rx')
for k=1:length(peakinds)
    text(y(peakinds(k),1),y(peakinds(k),2),num2str(k))
end
axis equal

subplot(2,1,2)
hold on
vmy=vecmag(y(:,3:4))
plot(t,vmy,'k')
for k=1:length(lumps)
    tau=(t'-lumps(k).C)/lumps(k).S+.5;
    tau=max(min(tau,1),0);
    kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/lumps(k).S;
    plot(t,kappa*norm(lumps(k).L))
end
for k=1:length(peakinds)
    text(t(peakinds(k)),vmy(peakinds(k)),num2str(k))
end

