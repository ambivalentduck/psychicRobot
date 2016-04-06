clc
clear all

load ../Data/curlkick/curlkick1Y.mat

xu=1.55;
nrows=4;
colors=linspecer(7);
lw=2;

f=find(([trials.targetcat]~=0)&([trials.disturbcat]));

inds=[22, 20, 38, 17, 11];
T=f(inds(5));

figure(1)
set(gcf,'units','centimeters','Position',[2,2,8.8,15])
clf

y=trials(T).y;
t=trials(T).ty;
inds=trials(T).i0:trials(T).if+20;
[~,peaks]=findpeaks(-vecmag(y(inds,3:4)));
inds=inds(1):inds(peaks(end)-5);
t=t-t(inds(1));

[lumps,resid]=findLumps(t',y(:,3:4),inds);

lowT=min([lumps.C]-[lumps.S]/2);
highT=max([lumps.C]+[lumps.S]/2);
inds=find((t>=lowT)&(t<=highT));
y=y(inds,:);
tshift=t(inds(1));
t=t(inds)-tshift;

[~,peakinds]=findpeaks(vecmag(y(:,3:4)));

[val,maxind]=max(vecmag(y(:,3:4)));

%% Plot x

subplot(nrows,1,1)
hold on
plot(y(:,1),y(:,2),'k','linewidth',lw)
%Choose that halfway through 2 is the alignment point.
[~,lumporder]=sort([lumps.C]);
lumpL=vertcat(lumps(lumporder).L);
lumpC=vertcat(lumps(lumporder).C);
lumpS=vertcat(lumps(lumporder).S);
cum=[0 0];
for k=1:length(lumps)
    [~,Ci]=min(abs(t-(lumpC(k)-tshift)));
    yinterp=interp1(t,y(:,1:2),t(Ci));
    lumpX0=yinterp;
    plot(lumpX0(1)+lumpL(k,1)/2*[-1 1],lumpX0(2)+lumpL(k,2)/2*[-1 1],'color',colors(k,:),'linewidth',lw)
    cum=cum+lumpL(k,:);
end
axis equal
plot(y(maxind,1),y(maxind,2),'rx','markersize',10)
text(y(maxind,1)-0.005,y(maxind,2)+0.04,'C')
plot([-.1 0],.38*[1 1],'k','linewidth',3)
text(-.05,.375,'10 centimeters','verticalalignment','top','horizontalalignment','center')
axis off
xl=xlim;
yl=ylim;
text(xl(1)-.02,mean(yl),'Intent','verticalalignment','bottom','horizontalalignment','center','rotation',90)
%Needs axis off and a 5cm mark

%% Plot v

subplot(nrows,1,2)
hold on
vmy=vecmag(y(:,3:4));
plot(t,vmy,'k','linewidth',lw)
for k=1:length(lumps)
    tau=(t'-lumpC(k)+tshift)/lumpS(k)+.5;
    tau=max(min(tau,1),0);
    kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/lumpS(k);
    plot(t,kappa*norm(lumpL(k,:)),'color',colors(k,:),'linewidth',lw)
end
plot(t(maxind),vmy(maxind),'rx','markersize',10)
text(t(maxind)-0.02,vmy(maxind)+.2,'C')
ylabel('Speed, meters/second')
set(gca,'ytick',[0 1])
xlim([0, xu])

%% Plot psi
subplot(nrows,1,3)
hold on
dr=y(peakinds(2),3:4);
dr=dr/norm(dr);
for kk=1:length(t)
    dp(kk)=dot(y(kk,3:4),dr)/norm(y(kk,3:4));
end
mingradpeakheight=.0001;
peak=peakinds(2);
gu=gradient(dp);
[~,locs1]=findpeaks(gu,'minpeakheight',mingradpeakheight);
[~,locs2]=findpeaks(-gu,'minpeakheight',mingradpeakheight);
lower=locs1(find(locs1<peak,1,'last'));
upper=locs2(find(locs2>peak,1,'first'));
plot(t,gu,'linewidth',lw,'color',colors(end,:))
plot(t(maxind),gu(maxind),'rx','markersize',10)
text(t(maxind)-0.02,gu(maxind)+.5,'C')
plot(t([lower maxind]),-.05+[0 0],'r')
text(mean(t([lower maxind])),-.05,'S/2','verticalalignment','top','horizontalalignment','center')
plot(t(lower),gu(lower),'ro','markersize',10)
ylabel('\psi','interpreter','tex')
xlim([0, xu])
xlabel('Time, seconds')

%% Plot x - y1

subplot(nrows,1,4)
hold on
tau=max(min((t-lumpC(3))/lumpS(3)+.5,1),0);
kappa=10*tau.^3-15*tau.^4+6*tau.^5;
yrec=y(:,1:2)-(lumpL(3,:)'*kappa)';


plot(yrec(:,1),yrec(:,2),'k','linewidth',lw)
%Choose that halfway through 2 is the alignment point.
[~,lumporder]=sort([lumps.C]);
lumpL=vertcat(lumps(lumporder).L);
lumpC=vertcat(lumps(lumporder).C);
lumpS=vertcat(lumps(lumporder).S);
cum=[0 0];
lumpC(2)=lumpC(2)+.325;
for k=1:length(lumps)
    if k==3
        continue
    end
    [~,Ci]=min(abs(t-(lumpC(k)-tshift)));
    yinterp=interp1(t,yrec(:,1:2),t(Ci));
    lumpX0=yinterp;
    plot(lumpX0(1)+lumpL(k,1)/2*[-1 1],lumpX0(2)+lumpL(k,2)/2*[-1 1],'color',colors(k,:),'linewidth',lw)
    cum=cum+lumpL(k,:);
end
axis equal
plot([-.1 0],.38*[1 1],'k','linewidth',3)
text(-.05,.375,'10 centimeters','verticalalignment','top','horizontalalignment','center')
%ylabel('Remainder')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
xl=xlim;
yl=ylim;
axis off
text(xl(1)-.02,mean(yl),'Remainder','verticalalignment','bottom','horizontalalignment','center','rotation',90)
