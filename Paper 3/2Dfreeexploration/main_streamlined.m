clc
clear all
close all

if ~exist('free_exp.mat','file')
    mat=load('h_1fch_ps1_JL.dat');
    
    tinds=70649:95876; %This particular subject
    
    t=0:.005:(mat(tinds(end),11)-mat(tinds(1),11));
    x=mat(tinds,[1 2]);
    
    filtn=64;
    filtType='loess';
    x=[smooth(x(:,1),filtn,filtType) smooth(x(:,2),filtn,filtType)];
    
    v=[gradient(x(:,1)) gradient(x(:,2))]/.005;
    a=[gradient(v(:,1)) gradient(v(:,2))]/.005;
    save('free_exp.mat','t','x','v','a');
else
    load('free_exp.mat')
end

%% Set constants
numbins=40;

%% Compute Ut, Tt
Ut=cumtrapz(t,dot(v',a')');
Tt=.5*(v(:,1).^2+v(:,2).^2);

%% Bin and count Ut levels, plot this as a Boltzmann distribution

[Utcounts,Utbins]=hist(Ut,numbins);
Utcounts=Utcounts/length(t);
logUtcounts=log(Utcounts)';
kTUt=[logUtcounts ones(size(logUtcounts))]\-Utbins';
R2Ut=cov(logUtcounts,-Utbins)/(std(logUtcounts)*std(Utbins));

figure(1)
clf
subplot(2,3,1)
hold on
bar(Utbins,Utcounts)
plot(Utbins,exp((-Utbins-kTUt(2))/kTUt(1)),'r','linewidth',3)
title(['Boltzmann Distribution Fit R^2 = ',num2str(R2Ut(1,2)^2)])
xlabel('Potential Energy/Constant mass, (m/s)^2')
ylabel(['Fraction of Time Points (N=',num2str(length(t)),')'])
legend('Data','Fit')

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

[Ttcounts,Ttbins]=hist(Tt,numbins);
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

[speedcounts,speedbins]=hist(speed,numbins);
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
