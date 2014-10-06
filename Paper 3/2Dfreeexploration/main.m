clc
clear all

mat=load('h_1fch_ps1_JL.dat');

tinds=70649:95876; %This particular subject

t=0:.005:(mat(tinds(end),11)-mat(tinds(1),11));
x=mat(tinds,[1 2]);

filtn=64;
filtType='loess';
x=[smooth(x(:,1),filtn,filtType) smooth(x(:,2),filtn,filtType)];

v=[gradient(x(:,1)) gradient(x(:,2))]/.005;
a=[gradient(v(:,1)) gradient(v(:,2))]/.005;

Ut=cumtrapz(t,dot(v',a')');

figure(1)
clf
plot(t,log(sort(Ut)))

lower=min(x);
upper=max(x);
range=upper-lower;

nbins=30;
binsize=range/nbins;

U=cell(nbins+1);

for k=1:length(t)
    inds=floor((x(k,:)-lower)./binsize)+1;
    U{inds(1),inds(2)}(end+1)=Ut(k);
end


Uxmean=zeros(size(U));
Uxstd=zeros(size(U));
for k=1:size(U,1)
    for kk=1:size(U,2)
        if isempty(U{k,kk})
            Uxmean(k,kk)=NaN;
        else
            Uxmean(k,kk)=mean(U{k,kk});
            Uxstd(k,kk)=std(U{k,kk})/sqrt(length(U{k,kk}));
        end
    end
end

figure(2)
clf
hold on
[X,Y]=meshgrid(lower(1):binsize(1):upper(1),lower(2):binsize(2):upper(2));
surf(X,Y,Uxmean)
mesh(X,Y,Uxmean-1.95*Uxstd,'FaceColor','None')
mesh(X,Y,Uxmean+1.95*Uxstd,'FaceColor','None')
xlabel('Robot X')
ylabel('Robot Y')
zlabel('Virtual Potential, U')