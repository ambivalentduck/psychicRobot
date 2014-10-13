function [R2b,R2p,kT,wattage]=boltzmannsubjectplots(t,x,v,a,fnumber)

%% Set constants
numbins=40;
lpctile=.1; %Throw away the craziest outliers
upctile=99.9;

%% Compute Ut, Tt
Ut=cumtrapz(t,dot(v',a')');
Ut=Ut-min(Ut); %Start as zero just for ease of reference

Tt=.5*(v(:,1).^2+v(:,2).^2);

%% Bin and count Ut levels, plot this as a Boltzmann distribution

lUt=prctile(Ut,lpctile);
uUt=prctile(Ut,upctile);
[Utcounts,Utbins]=hist(Ut,linspace(lUt,uUt,numbins));
Utcounts=Utcounts/length(t);
logUtcounts=log(Utcounts)';
kTUt=[logUtcounts ones(size(logUtcounts))]\-Utbins';
kT=kTUt(1);
R2Ut=cov(logUtcounts,-Utbins)/(std(logUtcounts)*std(Utbins));
R2b=R2Ut(1,2)^2;

figure(fnumber)
clf
subplot(2,3,1)
hold on
bar(Utbins,Utcounts)
plot(Utbins,exp((-Utbins-kTUt(2))/kTUt(1)),'r','linewidth',3)

power=abs(gradient(Ut,t));
lpower=prctile(power,lpctile);
upower=prctile(power,upctile);
[powercounts,powerbins]=hist(power,linspace(lpower,upower,numbins));
powercounts=powercounts/sum(powercounts);
logpowercounts=log(powercounts)';
recoeff=[logpowercounts ones(numbins,1)]\-powerbins';
wattage=recoeff(1);
R2Utmb=cov(logpowercounts,powerbins)/(std(logpowercounts)*std(powerbins));
R2p=R2Utmb(1,2)^2;
%plot(subUt,pmboltz,'g.')

title(['Boltzmann Distribution Fit R^2 = ',num2str(R2Ut(1,2)^2)])
xlabel('Potential Energy/Constant mass, (m/s)^2')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data','Boltzmann Fit','Maxwell-Boltzmann Fit')

%% Plot T vs U
subplot(2,3,2)
R2Lagrange=cov(Ut,Tt)/(std(Ut)*std(Tt));
plot(Ut,Tt,'k.','Markersize',.001)
mUt=max(Ut);
xlabel('Potential Energy/Constant mass, (m/s)^2')
ylabel('Kinetic Energy/Constant mass, (m/s)^2')
title(['Lagrangian Fit R^2 = ',num2str(R2Lagrange(1,2)^2)])
axis equal
xlim([0 mUt])
ylim([0 mUt])


%% This implies the relationship between v and P(v)

[theta,speed]=cart2pol(v(:,1),v(:,2));

[Ttcounts,Ttbins]=hist(Tt,linspace(prctile(Tt,lpctile),prctile(Tt,upctile),numbins));
Ttcounts=Ttcounts/length(t);
logTtcounts=log(Ttcounts)';
kTTt=[logTtcounts ones(size(logTtcounts))]\-Ttbins';
R2Tt=cov(logTtcounts,-Ttbins)/(std(logTtcounts)*std(Ttbins));

subplot(2,3,4)
hold on
bar(sqrt(2*Ttbins),Ttcounts)
plot(sqrt(2*Ttbins),exp((-Ttbins-kTTt(2))/kTTt(1)),'r','linewidth',3)
title(['Boltzmann Distribution Fit R^2 = ',num2str(R2Tt(1,2)^2)])
xlabel('Speed, m/s')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data','Fit')

[speedcounts,speedbins]=hist(speed,linspace(prctile(speed,lpctile),prctile(speed,upctile),numbins));
speedcounts=speedcounts/length(t);
logspeedcounts=log(speedcounts)';
kTspeed=[logspeedcounts ones(size(logspeedcounts))]\-(.5*speedbins.^2)';
R2speed=cov(logspeedcounts,-(.5*speedbins.^2))/(std(logspeedcounts)*std(-(.5*speedbins.^2)));

subplot(2,3,5)
hold on
bar(speedbins,speedcounts)
plot(speedbins,exp((-(.5*speedbins.^2)-kTspeed(2))/kTspeed(1)),'r','linewidth',3)
title(['Boltzmann Distribution Fit R^2 = ',num2str(R2speed(1,2)^2)])
xlabel('Speed, m/s')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data','Fit')

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

subplot(2,3,3)
contourf(thetas,speeds,flatcounts')
ylabel('Speed')
xlabel('Direction, radians')
title('Distribution of speed in direction')
colorbar

subplot(2,3,6)
[tgrid,sgrid]=meshgrid(thetas,speeds);
contourf(thetas,speeds,flatcounts'-exp((-(.5*sgrid.^2)-kTTt(2))/kTTt(1))/(numbins+1))
ylabel('Speed')
xlabel('Direction, radians')
title('Discrepancy in distribution of speed in direction')
colorbar
