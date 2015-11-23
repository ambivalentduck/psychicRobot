clc
clear all

empcolor='r';
fitcolor='b';
linewidth=1;

load free_exp_05stroke.mat

dt=mean(diff(t));
f=find((vecmag(v)>.25)|(vecmag(a)>1.5));
t=t(f);
x=x(f,:);
v=v(f,:);
speed=vecmag(v);
a=a(f,:);

%% Freeze time, vary distance

tspan=3; %seconds
indspan=round(tspan/dt);
nspans=floor(length(t)/indspan)

dists=zeros(nspans,1);

for k=1:nspans
    slide=1+(k-1)*indspan;
    dists(k)=sum(speed(slide:(slide+indspan))).^2;
end

figure(1)
clf
hold on

[Fd,xd,ld,ud]=ecdf(dists);
h=fill([xd(~isnan(ld)); wrev(xd(~isnan(ud)))],[ld(~isnan(ld)); wrev(ud(~isnan(ud)))],empcolor);
set(h,'edgealpha',0,'facealpha',.7)
plot(xd,Fd,'-','color',empcolor,'linewidth',linewidth)

[shift,n,A]=fitShiftedGam(dists)
plot(xd,gamcdf(xd-shift,n,A),'color',fitcolor,'linewidth',linewidth)


%% Freeze distance, vary time
