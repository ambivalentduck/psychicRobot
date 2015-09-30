clc
clear all

mass=.5; %.5 This affects nothing *but* plot magnitudes

if ~exist('cartoonFigData.mat','file')
    x1=[];
    while length(x1)~=4
        [x1,x2,t1,t2]=makeFreeSpaghetti;
    end
    save('cartoonFigData.mat','x1','x2','t1','t2')
else
    load('cartoonFigData.mat')
end

t1=t1/max(t2);
t2=t2/max(t2);

t=linspace(0,max(t2),200)';
dT=mean(diff(t));
vtotal=zeros(size(t));
atotal=zeros(size(t));

ptotal=zeros(size(t));
etotal=zeros(size(t));

for k=1:length(t1)
    inds=find((t>=t1(k))&(t<=t2(k)));
    tc=(t1(k)+t2(k))/2;
    ts=t2(k)-t1(k);
    ta=((t(inds)-tc))/ts+.5;
    kern=(x2(k)-x1(k))*(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
    acc=(x2(k)-x1(k))*(60*ta-180*ta.^2+120*ta.^3)/(ts^2);
    any(acc<0)
    
    subs(k).inds=inds;
    subs(k).t=t(inds);
    subs(k).v=kern;
    subs(k).a=acc;
    
    % Cartesian superposition
    vtotal(inds)=vtotal(inds)+kern;
    atotal(inds)=atotal(inds)+acc;
    
    % Energy superposition
    ptotal(inds)=ptotal(inds)+mass*kern.*acc;
    etotal(inds)=etotal(inds)+.5*mass*kern.^2;
end

cartP=mass*atotal.*vtotal;
cartE=.5*mass*vtotal.^2;
cartV=vtotal;

energyP=ptotal;
energyE=etotal;
energyV=sqrt(2*etotal/mass);

Nsub=3;
Msub=4;

cartColor='r';
energyColor='g';
subColor='k';
residSize=.4;
totalWidth=2;

figure(1)
clf

%% Power
subplot(Nsub,Msub,1)
hold on
plot(t,cartP,'-','color',cartColor,'linewidth',totalWidth)
plot(t,energyP,'-','color',energyColor,'linewidth',totalWidth)
resid=cartP;
for k=1:length(subs)
    plot(subs(k).t,mass*subs(k).a.*subs(k).v,'-','color',subColor)
    resid(subs(k).inds)=resid(subs(k).inds)-mass*subs(k).a.*subs(k).v;
end
plot(t,resid,'.','color',cartColor,'markersize',residSize)

subplot(Nsub,Msub,2)
tk=0:.001:1;
K=(60*tk-180*tk.^2+120*tk.^3).*(30*tk.^2-60*tk.^3+30*tk.^4);
[F,x]=ecdf(K);
plot(x,F,'k')

Nrand=10000;
rnd=rand(Nrand,2);
subplot(Nsub,Msub,4)
P=expinv(rnd(:,1),1).*(60*rnd(:,2)-180*rnd(:,2).^2+120*rnd(:,2).^3).*(30*rnd(:,2).^2-60*rnd(:,2).^3+30*rnd(:,2).^4);
ecdf(P)


%% Energy
subplot(Nsub,Msub,Msub+1)
hold on
plot(t,cartE,'-','color',cartColor,'linewidth',totalWidth)
plot(t,energyE,'-','color',energyColor,'linewidth',totalWidth)
resid=cartE;
for k=1:length(subs)
    plot(subs(k).t,.5*mass*subs(k).v.^2,'-','color',subColor)
    resid(subs(k).inds)=resid(subs(k).inds)-.5*mass*subs(k).v.^2;
end
plot(t,resid,'.','color',cartColor,'markersize',residSize)

subplot(Nsub,Msub,Msub+2)
Erange=0:.001:(.5*1.875^2);
PE=1-sqrt(1-2/15*sqrt(30)*sqrt(sqrt(2)*sqrt(Erange)));
plot(Erange,PE,'k')
xlim([-.1 .5*1.875^2+.1])

Nrand=10000;
rnd=rand(Nrand,2);
subplot(Nsub,Msub,Msub+4)
E=expinv(rnd(:,1),1).*(30*rnd(:,2).^2-60*rnd(:,2).^3+30*rnd(:,2).^4).^2;
ecdf(E)

%% Velocity
subplot(Nsub,Msub,2*Msub+1)
hold on
plot(t,cartV,'-','color',cartColor,'linewidth',totalWidth)
plot(t,energyV,'-','color',energyColor,'linewidth',totalWidth)
resid=energyV;
for k=1:length(subs)
    plot(subs(k).t,subs(k).v,'-','color',subColor)
    resid(subs(k).inds)=resid(subs(k).inds)-subs(k).v;
end
plot(t,abs(resid),'.','color',energyColor,'markersize',residSize)

subplot(Nsub,Msub,2*Msub+2)
vRange=0:.001:1.875;
Pv=1-sqrt(1-2*sqrt(30)*sqrt(vRange)/15);
plot(vRange,Pv,'k')
xlim([-.1 1.875+.1])

%% (L/T)^2
cdfarray=0:.1:5;
expcdfarray=expcdf(cdfarray,1);
for N=1:3
    subplot(Nsub,Msub,(N-1)*Msub+3)
    plot(cdfarray,expcdfarray,'k')
end

return
