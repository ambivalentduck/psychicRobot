clc
clear all

figure(1)
clf
set(gcf,'color','w')
Nsub=2;
Msub=2;

load free_exp_LP.mat

mass=.5; %.5 kg seems right, but don't know why.

%% Empirical CDFs
poscolor=[.2 .2 .8];
negcolor=[.8 .2 .2];
fitcolor=[.2 .8 .2];
linewidth=.5;

%Power
subplot(Nsub,Msub,1)
hold on
empiricalP=mass*dot(a',v')';

[F,p]=ecdf(abs(empiricalP));
%empiricalA=-1/(p(1:end-1)\log(-(F(1:end-1)-1)));
[empiricalA,ci]=expfit(abs(empiricalP))
pos=empiricalP(empiricalP>=0);
neg=-empiricalP(empiricalP<=0);
[Fp,xp,lp,up]=ecdf(pos);
[Fn,xn,ln,un]=ecdf(neg);

h=fill([p; wrev(p)],[expcdf(p,ci(1)); expcdf(wrev(p),ci(2))],fitcolor);
set(h,'edgealpha',0,'facealpha',.5)
h=fill([xp(~isnan(lp)); wrev(xp(~isnan(up)))],[lp(~isnan(lp)); wrev(up(~isnan(up)))],poscolor);
set(h,'edgealpha',0,'facealpha',.5)
h=fill([xn(~isnan(ln)); wrev(xn(~isnan(un)))],[ln(~isnan(ln)); wrev(un(~isnan(un)))],negcolor);
set(h,'edgealpha',0,'facealpha',.5)
ph=plot(xp,Fp,'-','color',poscolor,'linewidth',linewidth);
nh=plot(xn,Fn,'-','color',negcolor,'linewidth',linewidth);
fh=plot(p,expcdf(p,empiricalA),'-','color',fitcolor,'linewidth',linewidth);
xlim([-.01 2])
set(gca,'ytick',[0 1])
set(gca,'xtick',[0 1 2])
ylabel('CDF')
xlabel('Power, Watts')
legend([ph nh fh],{'+','-','Fit'},'Location','East')
%set(gca,'color','k')

