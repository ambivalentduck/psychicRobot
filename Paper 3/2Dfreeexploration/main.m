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
surf(Xgrid,Ygrid,Uxmean,Uxmean,'EdgeAlpha',0)
plot3(x(:,1),x(:,2),Ut,'k')
plot3(x(:,1),x(:,2),abs(gradient(Ut,t)),'r')
xlabel('Robot X')
ylabel('Robot Y')
zlabel('Virtual Potential, U')
grid on
set(gca,'CameraPosition',campos)

figure(2)
clf
hold on
[counts,bins]=hist(abs(gradient(Ut,t)),20);
counts=counts/sum(counts);
bar(bins,counts)

mb=[log(counts') ones(20,1)]\bins'

figure(99)
clf
hold on
plot(bins,counts,'.')
r2=cov(bins,log(counts))/(std(bins)*std(log(counts)))
plot(bins,exp((bins-mb(2))/mb(1)))
%plot(bins,exp((bins+mb(2))/mb(1)))

