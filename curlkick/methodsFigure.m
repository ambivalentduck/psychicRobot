clc
clear all

load ../Data/curlkick/curlkick1Y.mat

nrows=4;
colors=[1 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    0 1 1];

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

%% Plot x

subplot(nrows,1,1)
hold on
plot(y(:,1),y(:,2),'k')
%Choose that halfway through 2 is the alignment point.
[~,lumporder]=sort([lumps.C]);
lumpL=vertcat(lumps(lumporder).L);
lumpC=vertcat(lumps(lumporder).C);
lump0i=[4];
cum=[0 0];
lumpX0=[0 0];
for k=1:length(lumps)
    if any(k==lump0i)
        [~,Ci]=min(abs(t-(lumpC(k)-tshift)));
        yinterp=interp1(t,y(:,1:2),t(Ci));
        plot(y(Ci,1),y(Ci,2),'mo')
        (cum+lumpL(k,:)/2-yinterp)
        lumpX0=lumpX0+(cum+lumpL(k,:)/2-yinterp)/length(lump0i);
    end
    cum=cum+lumpL(k,:);
end
%lumpX0=-(lumpL(1,:)+lumpL(2,:)/2-y(peakinds(2),1:2));
lumpVertices=cumsum([-lumpX0;lumpL]);
plot(lumpVertices(:,1),lumpVertices(:,2),'m')

%plot(y(peakinds,1),y(peakinds,2),'rx')

for k=1:length(peakinds)
    text(y(peakinds(k),1),y(peakinds(k),2),['C_',num2str(k)])
end
axis equal

%Needs axis off and a 5cm mark

%% Plot v

subplot(nrows,1,2)
hold on
vmy=vecmag(y(:,3:4));
plot(t,vmy,'k')
for k=1:length(lumps)
    tau=(t'-lumps(k).C+tshift)/lumps(k).S+.5;
    tau=max(min(tau,1),0);
    kappa=(30*tau.^2-60*tau.^3+30*tau.^4)/lumps(k).S;
    plot(t,kappa*norm(lumps(k).L))
end
for k=1:length(peakinds)
    text(t(peakinds(k)),vmy(peakinds(k)),num2str(k))
end

subplot(nrows,1,3)
dr=y(peakinds(2),3:4);
dr=dr/norm(dr);
for kk=1:length(t)
    dp(kk)=dot(y(kk,3:4),dr)/norm(y(kk,3:4));
end
plot(t,dp)
ylabel('$\hat{L}_1 \cdot \hat{v}$','interpreter','latex')

subplot(nrows,1,4)
mingradpeakheight=.0001;
peak=peakinds(2);
    gu=gradient(dp);
    [~,locs1]=findpeaks(gu,'minpeakheight',mingradpeakheight);
    [~,locs2]=findpeaks(-gu,'minpeakheight',mingradpeakheight);
    lower=locs1(find(locs1<peak,1,'last'));
    upper=locs2(find(locs2>peak,1,'first'));
    plot(t,gu)


