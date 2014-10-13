clc
clear all
close all

campos=[3.4800   -3.7304    0.0445];

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

%% Compute Ut
Ut=cumtrapz(t,dot(v',a')');

%% Compute Ux
lower=min(x);
upper=max(x);
range=upper-lower;

nbins=30;
binsize=range/nbins;

[Xgrid,Ygrid]=meshgrid(lower(1):binsize(1):upper(1),lower(2):binsize(2):upper(2));

U=cell(nbins+1);
for k=1:length(t)
    inds=floor((x(k,:)-lower)./binsize)+1;
    U{inds(1),inds(2)}(end+1)=Ut(k);
end

Uxmean=zeros(size(U));
Uxstd=zeros(size(U));
counts=cell(size(U));
flatcounts=zeros(size(U));
for k=1:size(U,1)
    for kk=1:size(U,2)
        if isempty(U{k,kk})
            Uxmean(k,kk)=NaN;
        else
            Uxmean(k,kk)=mean(U{k,kk});
            l=length(U{k,kk});
            Uxstd(k,kk)=std(U{k,kk})/sqrt(l);
            counts{k,kk}=l*ones(l,1);
            flatcounts(k,kk)=l;
        end
    end
end


%% Show how we’re doing this and how consistent it is or isn't.
FIG=0;
FIG=FIG+1;
figure(FIG)
clf
hold on
surf(Xgrid,Ygrid,Uxmean,log(Uxmean),'EdgeAlpha',0)
plot3(x(:,1),x(:,2),Ut,'k')
xlabel('Robot X')
ylabel('Robot Y')
zlabel('Virtual Potential, U')
grid on
set(gca,'CameraPosition',campos)

FIG=FIG+1;
figure(FIG)
clf
hold on
plot3(v(:,1),v(:,2),Ut,'k')
xlabel('Velocity, Robot X/s')
ylabel('Velocity, Robot Y/s')
zlabel('Virtual Potential, U')
grid on
set(gca,'CameraPosition',campos)

FIG=FIG+1;
figure(FIG)
clf
kin=.5*vecmag(v).^2;
R2=cov(kin,Ut)/(std(kin)*std(Ut))
plot(kin,Ut,'.','Markersize',.001)
xlabel('Speed^2: Kinetic Energy, (m/s)^2')
ylabel('U, Energy Units')
title('Kinetic vs Potential Energy. Both are scaled by the same mass.')

%%


%% Each bin in Ux “owns” a count and a mean U, this is 1D and lets us fit T
FIG=FIG+1;
figure(FIG)
clf
%subplot(1,2,1) %Whole data, unfair since what's P[x] really?
hold on
C1D=vertcat(counts{:});
C1Dnorm=C1D/sum(C1D);
f=find(C1Dnorm>0);
lC1Dnorm=log(C1Dnorm(f));
U1D=horzcat(U{:});
U1D=U1D(f);
kToff=[lC1Dnorm,ones(size(lC1Dnorm))]\(-U1D');
plot(lC1Dnorm,-U1D,'.','Markersize',.001)
fitX=linspace(min(lC1Dnorm),max(lC1Dnorm),20);
plot(fitX,kToff(1)*fitX+kToff(2),'r')
R2=cov(lC1Dnorm,-U1D)/(std(lC1Dnorm)*std(-U1D))
title(['R^2 = ',num2str(R2(1,2)^2)])
ylabel('Calculated Potential, Energy Units')
xlabel('log Fraction of Time Points in a Potential Range')

% subplot(1,2,2) %Means and bin frequences is probably more informative
% hold on
% C1D=flatcounts(:);
% C1Dnorm=C1D/sum(C1D);
% f=find(C1Dnorm>0);
% lC1Dnorm=log(C1Dnorm(f));
% U1D=Uxmean(:);
% U1D=U1D(f);
% kToff=[lC1Dnorm,ones(size(lC1Dnorm))]\(-U1D);
% plot(lC1Dnorm,-U1D,'.','Markersize',.001)
% fitX=linspace(min(lC1Dnorm),max(lC1Dnorm),20);
% plot(fitX,kToff(1)*fitX+kToff(2),'r')
% R2=cov(lC1Dnorm,-U1D)/(std(lC1Dnorm)*std(-U1D))
% title(['R^2 = ',num2str(R2(1,2)^2)])
% ylabel('Calculated Potential, Energy Units')
% xlabel('log Fraction of Time Points in a Potential Range')

%% We have U and T, so now predict U in 2D, but the statistical test was already done above when we got our R^2
FIG=FIG+1;
figure(FIG)
clf
subplot(1,2,1)
hold on
flatcountsnorm=flatcounts/sum(sum(flatcounts));
surf(Xgrid,Ygrid,flatcountsnorm,'EdgeAlpha',0)
predictedcounts=exp((-Uxmean-kToff(2))/kToff(1));
predictedt=exp((-Ut-kToff(2))/kToff(1));
%predictedt=predictedt/sum(predictedt);
plot3(x(:,1),x(:,2),exp((-Ut-kToff(2))/kToff(1)),'k.','Markersize',.001)
grid on
set(gca,'CameraPosition',campos)

subplot(1,2,2)
hold on
surf(Xgrid,Ygrid,predictedcounts)
grid on
set(gca,'CameraPosition',campos)

%% There's structure, clearly, even if it's not position

FIG=FIG+1;
figure(FIG)
clf
hold on
[N,cent]=hist(Ut,20);
N=N/sum(N);
bar(cent,N)
kToff=[log(N)' ones(size(N))']\-cent';
plot(cent,exp((-cent-kToff(2))/kToff(1)),'r','linewidth',3)
R2=cov(log(N),-cent)/(std(log(N))*std(cent));
title(['R^2 = ',num2str(R2(1,2)^2)])
xlabel('Calculated Potential, Energy Units')
ylabel('Fraction of Time Points in a Potential Range')

%% Based on the lagrangian, all of the structure is in velocity, so let's prove it

lower=min(v);
upper=max(v);
range=upper-lower;

nbins=30;
binsize=range/nbins;

[Vxgrid,Vygrid]=meshgrid(lower(1):binsize(1):upper(1),lower(2):binsize(2):upper(2));

flatcounts=zeros(nbins+1);
for k=1:length(t)
    inds=floor((v(k,:)-lower)./binsize)+1;
    flatcounts(inds(1),inds(2))=flatcounts(inds(1),inds(2))+1;
end
flatcounts=flatcounts/length(t);

[N,cent]=hist(Ut,20);
N=N/sum(N);
bar(cent,N)
kToff=[log(N)' ones(size(N))']\-cent';

FIG=FIG+1;
figure(FIG)
clf
hold on
surf(Vxgrid,Vygrid,PU,'EdgeColor','k','FaceAlpha',0)
T=.5*(Vxgrid.^2+Vygrid.^2);
PU=exp((-T-kToff(2))/kToff(1))
PU=PU/sum(PU(:));
surf(Vxgrid,Vygrid,flatcounts,'EdgeAlpha',0,'FaceAlpha',.7)
xlabel('Robot Velocity X-component, m/s')
ylabel('Robot Velocity Y-component, m/s')
zlabel('Probability')

FIG=FIG+1;
figure(FIG)
clf
hold on
surf(Vxgrid,Vygrid,flatcounts-PU)
xlabel('Robot Velocity X-component, m/s')
ylabel('Robot Velocity Y-component, m/s')
zlabel('Difference in Probability')

%% So, re-represent this whole thing and graph it a bit differently
[theta,speed]=cart2pol(v(:,1),v(:,2));

FIG=FIG+1;
figure(FIG)
clf
subplot(1,2,1)
hold on

[N,cent]=hist(.5*speed.^2,20);
bar(sqrt(cent),N)

kToff=[log(N)' ones(size(N))']\-cent';
plot(sqrt(cent),exp((-cent-kToff(2))/kToff(1)),'r','linewidth',3)
R2=cov(log(N),-cent)/(std(log(N))*std(cent));
title(['R^2 = ',num2str(R2(1,2)^2)])
xlabel('Speed, m/s')
ylabel('Fraction of Time Points in a Speed Range')
legend('Actual','Fit')

subplot(1,2,2)


