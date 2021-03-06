clc
clear all

font='Arial';
fontsize=12;
axlabsize=10;
interpreter='tex';

figure(1)
clf
set(gcf,'color','w')
Nsub=3;
Msub=2;

spHandles=zeros(Nsub*Msub,1);

load free_exp_LP.mat

mass=.5; %.5 kg seems right, but don't know why.

%% Empirical CDFs
poscolor=[.2 .2 .8];
negcolor=[.8 .2 .2];
fitcolor=[.2 .8 .2];
linewidth=.5;

%Power
spHandles(1)=subplot(Nsub,Msub,1);
hold on
empiricalP=mass*dot(a',v')';

[F,p]=ecdf(abs(empiricalP));
%empiricalA=-1/(p(1:end-1)\log(-(F(1:end-1)-1)));
[empiricalA,ci]=expfit(abs(empiricalP));
pos=empiricalP(empiricalP>=0);
neg=-empiricalP(empiricalP<=0);
[Fp,xp,lp,up]=ecdf(pos);
[Fn,xn,ln,un]=ecdf(neg);

h=fill([p; wrev(p)],[expcdf(p,ci(1)); expcdf(wrev(p),ci(2))],fitcolor);
set(h,'edgealpha',0,'facealpha',.5);
h=fill([xp(~isnan(lp)); wrev(xp(~isnan(up)))],[lp(~isnan(lp)); wrev(up(~isnan(up)))],poscolor);
set(h,'edgealpha',0,'facealpha',.5)
h=fill([xn(~isnan(ln)); wrev(xn(~isnan(un)))],[ln(~isnan(ln)); wrev(un(~isnan(un)))],negcolor);
set(h,'edgealpha',0,'facealpha',.5)
ph=plot(xp,Fp,'-','color',poscolor,'linewidth',linewidth);
nh=plot(xn,Fn,'-','color',negcolor,'linewidth',linewidth);
fh=plot(p,expcdf(p,empiricalA),'-','color',fitcolor,'linewidth',linewidth);
xlim([-.05 2])
set(gca,'ytick',[0 1])
set(gca,'xtick',[0 2])
ylabel('CDF')
xlabel('Power, Watts')
legend([ph nh fh],{'+','-','Fit'},'Location','East')

%Energy
spHandles(3)=subplot(Nsub,Msub,3);
hold on
empiricalE=mass*dot(v',v');
[Fe,xe,le,ue]=ecdf(empiricalE);
h=fill([xe(~isnan(le)); wrev(xe(~isnan(ue)))],[le(~isnan(le)); wrev(ue(~isnan(ue)))],poscolor);
set(h,'edgealpha',0,'facealpha',.5)
plot(xe,Fe,'-','color',poscolor,'linewidth',linewidth)
xlim([-.01 .5])
set(gca,'ytick',[0 1])
set(gca,'xtick',[0 .5])
ylabel('CDF')
xlabel('Energy, Joules')
set(gca,'linewidth',2,'ticklength',[0;0])
%nA=gamfit(empiricalE);
%n=nA(1);
%A=nA(2);
%shift=0;
[shift,n,A]=fitShiftedGam(xe);
plot(xe,gamcdf(xe-shift,n,A),'color',fitcolor)

%% Cartoon hypotheses

x1=[0 .02 .08]*.50/.1181;
x2=[.03 .10 .10]*.50/.1181;
durs=[.5 .7 .6]/1.2;
t1=[0 0 0 0];
t2=[0 0 0];
for k=1:length(durs)
    t2(k)=t1(k)+durs(k);
    t1(k+1)=t1(k)+.5*durs(k);
end
t1(end)=[];

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
    
    subs(k).inds=inds; %#ok<*SAGROW>
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

cartColor='r';
energyColor='g';
subColor='k';
residSize=1;
totalWidth=2;

%Power
spHandles(2)=subplot(Nsub,Msub,2);
hold on
hcart=plot(t,cartP,'-','color',cartColor,'linewidth',totalWidth);
henergy=plot(t,energyP,'-','color',energyColor,'linewidth',totalWidth);
resid=cartP;
for k=1:length(subs)
    hlump=plot(subs(k).t,mass*subs(k).a.*subs(k).v,'-','color',subColor);
    resid(subs(k).inds)=resid(subs(k).inds)-mass*subs(k).a.*subs(k).v;
end
plot(t,resid,'.','color',cartColor,'markersize',residSize)
ylim([-2 2])
set(gca,'ytick',[-2 0 2])
set(gca,'xtick',[0 1])
ylabel('Power, Watts')
text(.13,.8,'\{','rotation',-90,'horizontalalignment','center','verticalalignment','middle','fontsize',24)
text(.128,1.15,'1','horizontalalignment','center','fontsize',12)
text(.27,-.9,'\{','rotation',90,'horizontalalignment','center','verticalalignment','middle','fontsize',24)
text(.28,-1.3,'2','horizontalalignment','center','fontsize',12)
legend([hlump,henergy,hcart],{'Subunit','\Sigma Energy','\Sigma Velocity'})

%Energy
spHandles(4)=subplot(Nsub,Msub,4);
hold on
plot(t,cartE,'-','color',cartColor,'linewidth',totalWidth)
plot(t,energyE,'-','color',energyColor,'linewidth',totalWidth)
resid=cartE;
for k=1:length(subs)
    plot(subs(k).t,.5*mass*subs(k).v.^2,'-','color',subColor)
    resid(subs(k).inds)=resid(subs(k).inds)-.5*mass*subs(k).v.^2;
end
plot(t,resid,'.','color',cartColor,'markersize',residSize)
set(gca,'ytick',[0 .3])
ylim([0 .3])
set(gca,'xtick',[0 1])
ylabel('Energy, Joules')
xlabel('Time, seconds')

% Velocity
spHandles(6)=subplot(Nsub,Msub,6)
hold on
plot(t,cartV,'-','color',cartColor,'linewidth',totalWidth)
plot(t,energyV,'-','color',energyColor,'linewidth',totalWidth)
resid=energyV;
for k=1:length(subs)
    plot(subs(k).t,subs(k).v,'-','color',subColor)
    resid(subs(k).inds)=resid(subs(k).inds)-subs(k).v;
end
plot(t,abs(resid),'.','color',energyColor,'markersize',residSize)
ylabel('Velocity, meters/second')
xlabel('Time, seconds')

%% Quick sim
N=floor(10*60/.3);
L2Tn2=exprnd(.09,N,1);
Tn2=1*ones(N,1); %4.6+exprnd(7,N,1);
dx=sqrt(L2Tn2./Tn2);
t1=zeros(N,1);
t2=t1;
for k=1:N
   t2(k)=t1(k)+1/sqrt(Tn2(k));
   xp=.3; %This works
   if xp>1
       xp=1;
   end
   if k~=N
       t1(k+1)=t1(k)+xp*(t2(k)-t1(k));
   end
end
t=0:.005:t2(end);
vtotal=zeros(size(t));
atotal=zeros(size(t));
etotal=zeros(size(t));
ptotal=zeros(size(t));

for k=1:N
    inds=find((t>=t1(k))&(t<=t2(k)));
    tc=(t1(k)+t2(k))/2;
    ts=t2(k)-t1(k);
    ta=((t(inds)-tc))/ts+.5;
    kern=dx(k)*(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
    acc=dx(k)*(60*ta-180*ta.^2+120*ta.^3)/(ts^2);
    
    % Cartesian superposition
    vtotal(inds)=vtotal(inds)+kern;
    atotal(inds)=atotal(inds)+acc;
    
    % Energy superposition
    ptotal(inds)=ptotal(inds)+mass*kern.*acc;
    etotal(inds)=etotal(inds)+.5*mass*kern.^2;
end
cartPm=mass*atotal.*vtotal;
cartEm=.5*mass*vtotal.^2;
cartVm=vtotal;

energyPm=ptotal;
energyEm=etotal;
energyVm=sqrt(2*etotal/mass);

spHandles(5)=subplot(Nsub,Msub,5);
hold on
empiricalE=mass*dot(v',v');
%[Fe,xe,le,ue]=ecdf(empiricalE);
[Fcm,xcm,lcm,ucm]=ecdf(cartEm);
[Fem,xem,lem,uem]=ecdf(energyEm);
plot(xe,Fe,'b')
plot(xem,Fem,'g')
plot(xcm,Fcm,'r')
xlim([-.1 .5])

expfit(abs(empiricalP))
gamfit(empiricalE)
expfit(abs(energyPm))
gamfit(energyEm)

%% Cleanup at the end

%Setup axes
allaxes=findobj(gcf,'type','axes');
set(allaxes,'linewidth',2,'ticklength',[0;0]);
set(gcf,'units','inches')
figW=8.5;
figH=11/2;
set(gcf,'position',[8 5 figW figH]);

lmargin=.08;
rmargin=.02;
tmargin=.02;
bmargin=.02;
rowHeights=[.25 .25 .25];
rowMargins=[0 .08 .08];
rowBottoms=1-cumsum(rowHeights+rowMargins)-tmargin;
colWidth=.38;
colX=[lmargin 1-rmargin-colWidth];
for k=1:length(spHandles)
    if spHandles(k)==0
        continue
    end
    column=2-mod(k,2);
    row=floor((k-1)/2)+1;
    set(spHandles(k),'position',[colX(column) rowBottoms(row) colWidth rowHeights(row)]);
end

%Format all text
alltext=findobj(gcf,'type','text');

textInvariant=['\fontname{',font,'} \fontsize{',num2str(fontsize),'} '];
for k=1:length(alltext)
    string=get(alltext(k),'string');
    set(alltext(k),'string',[textInvariant, string],'interpreter',interpreter);
end


