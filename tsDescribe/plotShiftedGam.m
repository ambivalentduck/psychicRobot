function plotShiftedGam(Y)

hold on

[F,x,l,u]=ecdf(Y);
h=fill([x(~isnan(l)); wrev(x(~isnan(u)))],[l(~isnan(l)); wrev(u(~isnan(u)))],'k');
set(h,'edgealpha',0,'facealpha',.5);
plot(x,F,'k','linewidth',2)

[shift,n,A]=fitShiftedGam(Y);
plot(x,gamcdf(x-shift,n,A),'color','r','linewidth',2)

mu=mean(Y);
sigma=std(Y);
plot(x,normcdf(x,mu,sigma),'b','linewidth',2)
