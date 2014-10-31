function [R2U,R2P,T,W]=boltzmannsubjectplots(t,x,v,a,fnumber)

global fitmecounts fitmebins

%% Set constants
numbins=40;
lpctile=10;
upctile=99.9; %Throw away the craziest outliers to help hist not find zeros

sWidth=2;
sHeight=4;

%% Compute Ut, Tt
Pt=dot(v',a')';
Ut=cumtrapz(t,Pt);
Ut=Ut-min(Ut); %Flush to zero just for ease of reference
Tt=.5*(v(:,1).^2+v(:,2).^2);

%% Bin and count Ut levels, plot this as a Boltzmann distribution

figure(fnumber)
clf
subplot(sWidth,sHeight,1)
[R2P,W,aPbins,aPcounts,Pfitcounts]=histBoltzmann(abs(Pt),lpctile,upctile,1);
title(['Boltzmann Distribution Fit R^2 = ',num2str(R2P,3)])
xlabel('|Work|/Mass, (m/s)^2')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data',['Boltzmann Fit, W=',num2str(W,2)])

subplot(sWidth,sHeight,2)
hold on
[Pcounts,Pbins]=hist(Pt,64);
Pcounts=Pcounts/length(t);
bar(Pbins,Pcounts)
plot(Pbins,1/(2*W)*exp(-abs(Pbins)/W),'r','linewidth',2)

title(['Boltzmann Distribution Fit R^2 = ',num2str(R2P,3)])
xlabel('Work/Mass, (m/s)^2')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data',['Boltzmann Fit, W=',num2str(W,2)])

subplot(sWidth,sHeight,3)
hold on
[R2U,T,Ubins,Ucounts,Ufitcounts]=histBoltzmann(Ut,lpctile,upctile,1);
nmin=find(Ubins>0,1,'first');
Usubbins=Ubins(nmin:end);
Usub=Ut((Ut>=Usubbins(1))&(Ut<=Usubbins(end)));
Usubcounts=hist(Usub,Usubbins);
Usubcounts=Usubcounts/sum(Usubcounts);
fitmecounts=Usubcounts;
fitmebins=Usubbins;

xfit=fminunc(@bruteforce,[2.3 .2])

gam=gampdf(Ubins,xfit(1),xfit(2));
gam=gam/sum(gam);
plot(Ubins,gam,'g','linewidth',3)
1-var(Ucounts-gam)/var(Ucounts)
ylim([0 1.2*max(Ucounts)])
title(['Boltzmann Distribution Fit R^2 = ',num2str(R2U,3)])
xlabel('Potential Energy/Mass, (m/s)^2')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data',['Boltzmann Fit, T=',num2str(T,2)],['Gamma Dist, n=',num2str(xfit(1),3),', T=',num2str(xfit(2),2)])


%% Plot T vs U
% subplot(sWidth,sHeight,5)
% R2Lagrange=cov(Ut,Tt)/(std(Ut)*std(Tt));
% plot(Ut,Tt,'k.','Markersize',.001)
% mUt=max(Ut);
% xlabel('Potential Energy/Mass, (m/s)^2')
% ylabel('Kinetic Energy/Mass, (m/s)^2')
% title(['Lagrangian Fit R^2 = ',num2str(R2Lagrange(1,2)^2)])
% axis equal
% xlim([0 mUt])
% ylim([0 mUt])

%% Fit the relationship between v and P(v)

[theta,speed]=cart2pol(v(:,1),v(:,2));

upper=prctile(speed,upctile);
lower=prctile(speed,lpctile);

[speedcounts,speedbins]=hist(speed((speed<=upper)&(speed>=lower)),linspace(lower,upper,64));
speedcounts=speedcounts/sum(speedcounts);

subplot(sWidth,sHeight,5)
hold on
bar(speedbins,speedcounts)

fitmecounts=speedcounts;
fitmebins=.5*speedbins.^2;
xfit=fminunc(@bruteforce,[2.5 .2])
xfitspeed=xfit;
gam=gampdf(.5*speedbins.^2,xfit(1),xfit(2));
gam=gam/sum(gam);
plot(speedbins,gam,'g','linewidth',3)

xlabel('Speed, m/s')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data',['Gamma Dist, n=',num2str(xfit(1),3),', T=',num2str(xfit(2),2)])

%% Speed and direction are actually orthogonal in polar coordinates

polcord=[theta,speed];
lower=min(polcord);
upper=max(polcord);
binsize=(upper-lower)./(numbins*[1 1]);

flatcounts=zeros(numbins+1);
for k=1:length(t)
    inds=floor((polcord(k,:)-lower)./binsize)+1;
    flatcounts(inds(1),inds(2))=flatcounts(inds(1),inds(2))+1;
end
flatcounts=flatcounts/length(t);

speeds=linspace(lower(2),upper(2),numbins+1);
thetas=linspace(lower(1),upper(1),numbins+1);

subplot(sWidth,sHeight,4)
contourf(thetas,speeds,flatcounts')
ylabel('Speed')
xlabel('Direction, radians')
title('Distribution of speed in direction')
colorbar

subplot(sWidth,sHeight,8)
[tgrid,sgrid]=meshgrid(thetas,speeds);
gam=gampdf(.5*sgrid.^2,xfit(1),xfit(2));
gam=gam/sum(gam(:));
contourf(thetas,speeds,flatcounts'-gam)
ylabel('Speed')
xlabel('Direction, radians')
title('Discrepancy in distribution of speed in direction')
colorbar

subplot(sWidth,sHeight,6:7)
mins=min(x);
spans=max(x)-mins;
nbins=64;
xhist=zeros(nbins);
a1hist=zeros(nbins);
a2hist=zeros(nbins);

for k=1:length(t)
    inds=floor(nbins*(x(k,:)-mins)./spans);
    inds=min(inds+1,[nbins nbins]);
    xhist(inds(1),inds(2))=xhist(inds(1),inds(2))+1;
    a1hist(inds(1),inds(2))=a1hist(inds(1),inds(2))+a(k,1);
    a2hist(inds(1),inds(2))=a2hist(inds(1),inds(2))+a(k,2);
end
a1hist=a1hist./xhist;
a2hist=a2hist./xhist;
xhist=xhist/length(t);

smoothxhist=filter2(ones(6)/36,xhist);

U=-T*log(T*smoothxhist);
U(U==inf)=-pi;
U(U==-pi)=max(U(:));
%[a1hist,a2hist]=gradient(U);


[gradx1,gradx2]=gradient(smoothxhist);
maggradx=sqrt(gradx1.^2+gradx2.^2);
maggradU=sqrt(a1hist.^2+a2hist.^2);
riseorfall=(a1hist.*gradx1+a2hist.*gradx2)./(maggradx.*maggradU);
riseorfall(isnan(riseorfall))=-1;

gamU=linspace(0,max(Ut),nbins);
gam=gampdf(gamU,xfitspeed(1),xfitspeed(2));
[maxgam,imax]=max(gam)
xhist=maxgam*xhist/max(xhist(:));
max(xhist(:))

risingP=gam(1:imax);
risingU=gamU(1:imax);

fallingP=gam(imax:end);
fallingU=gamU(imax:end);

risingVirtualU=interp1(risingP,risingU,xhist);
fallingVirtualU=interp1([fallingP 0],[fallingU fallingU(end)*.9],xhist);

figure(fnumber+20)
[x1grid,x2grid]=meshgrid(linspace(mins(1),mins(1)+spans(1),nbins),linspace(mins(2),mins(2)+spans(2),nbins));
contourf(x1grid,x2grid,U)
xlabel('Robot X-axis, m')
ylabel('Robot Y-axis, m')
title('Virtual Potential via Boltzmann')
axis equal
colorbar

figure(fnumber+21)
surf(x1grid,x2grid,gradient(smoothxhist))

figure(fnumber+22)
contourf(x1grid,x2grid,riseorfall)
colorbar

figure(fnumber+23)
contourf(x1grid,x2grid,risingVirtualU)
axis equal
colorbar
caxis([0 gamU(end)])

figure(fnumber+24)
contourf(x1grid,x2grid,fallingVirtualU)
axis equal
colorbar
caxis([0 gamU(end)])

figure(fnumber+25)
combined=(riseorfall>0).*risingVirtualU+(riseorfall<=0).*fallingVirtualU;
contourf(x1grid,x2grid,combined)
axis equal
colorbar
caxis([0 gamU(end)])

figure(fnumber+26)
contourf(x1grid,x2grid,U)
axis equal
colorbar
%caxis([0 gamU(end)])

error




function cost=bruteforce(x)
global fitmecounts fitmebins

gam=gampdf(fitmebins,x(1),x(2));
gam=gam/sum(gam);
cost=sum((fitmecounts-gam).^2);
