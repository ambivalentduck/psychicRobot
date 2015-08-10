clc
clear all
close all

load('free_exp_05stroke.mat')
%load('free_exp_MF.mat')

en=dot(v',v');
figure(74)
clf
inds=1:20:size(x,1);
plot3(x(inds,1),x(inds,2),en(inds),'.')

mu=expfit(en);

figure(75)
clf
<<<<<<< HEAD
sx1=size(x,1);
ds=1:150:sx1;
plot3(x(ds,1),x(ds,2),en(ds),'.')
=======

lower=min(x);
upper=max(x);
range=upper-lower;
nbins=64;
selectors=[(x(:,1)-lower(1))/range(1), (x(:,2)-lower(2))/range(2)];
selectors=floor(selectors*64)+1;
for k=1:2
    selectors(selectors(:,k)==65,k)=64;
end

counts=zeros(nbins);
for k=1:nbins
    f=find(selectors(:,1)==k);
    for kk=1:nbins
        counts(k,kk)=sum(selectors(f,2)==kk);
    end
    counts(k,counts(k,:)==0)=1;
end

U=-log(counts)*mu;
Unorm=U-min(U(:));
xbins=linspace(lower(1),upper(1),nbins);
ybins=linspace(lower(2),upper(2),nbins);
imagesc(xbins,ybins,Unorm)
axis equal
hold on
plot(-.25+[0 .05],-.6*[1 1],'k-')
text(-.25+.025, -.6,'5 cm','HorizontalAlignment','center','verticalalignment','top')
axis off


figure(76)
clf
hold on
[xgrid,ygrid]=meshgrid(xbins,ybins);
mesh(xgrid',ygrid',Unorm,'facealpha',0)
inds=1:15:size(x,1);
enNorm=en;
for k=1:length(en)
    enNorm(k)=en(k)+Unorm(selectors(k,1),selectors(k,2));
end
plot3(x(inds,1),x(inds,2),enNorm(inds),'.','markersize',10)
zlabel('Total Energy, J')
ylabel('Robot y, m')
xlabel('Robot x, m')

figure(77)
clf
hold on
[sun,ord]=sort(Unorm(:));
for k=1:length(ord)
        XR=mod(ord(k)-1,nbins)+1;
        YR=floor(ord(k)/nbins)+1;
        YR=min(64,YR);
        f=find((selectors(:,1)==XR)&(selectors(:,2)==YR));
        plot((k/(nbins^2))*ones(length(f),1),en(f)+Unorm(XR,YR),'.','markersize',.001)
end
plot(linspace(0,1,nbins^2),sun,'r-')
set(gca,'xtick',[])
xlabel('Potential Rank')
ylabel('Total Energy, J')
>>>>>>> origin/master
