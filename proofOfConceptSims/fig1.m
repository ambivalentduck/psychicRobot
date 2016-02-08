clc
clear all

dU=1;
L2=1;
Rt=1;

Tn2crit=dU/(3*L2*Rt);
Tcrit=Tn2crit^-.5;
Tn2_0=dU/(L2*Rt);

Sn2=linspace(0,dU/(L2*Rt),50);

dJ=dU-L2*Sn2*Rt;

figure(1)
clf
subplot(2,2,1)
plot(Sn2,dJ,'-')
ylabel('J')
xlim([0 1.05*Tn2_0])
set(gca,'xtick',[])

subplot(2,2,3)
S=Sn2.^-.5;
dJdot=dU./S-L2*Rt./(S.^3);
plot(Sn2,dJdot,'--.')
text(Tn2crit,dU/Tcrit-L2*Rt/(Tcrit^3),'$\frac{\Delta U}{3L^2R_T}$','interpreter','latex','fontsize',20,'horizontalalignment','center','verticalalignment','top')
text(1,-0.01,'$\frac{\Delta U}{L^2 R_T}$','interpreter','latex','fontsize',20,'horizontalalignment','center','verticalalignment','top')
xlabel('S^{-2}')
ylabel('J dot')
xlim([0 1.05*Tn2_0])
set(gca,'xtick',[0 Tn2_0])
set(gca,'xticklabels',{'0',''})

subplot(2,2,2)
plot(Sn2.^-.5,dJ,'-')
ylabel('J')
%xlim([0 1.05*Tn2_0])
set(gca,'xtick',[])

subplot(2,2,4)
plot(Sn2.^-.5,dJdot,'--.')
text(1,-0.01,'$\sqrt{\frac{L^2 R_T}{\Delta U}}$','interpreter','latex','fontsize',16,'horizontalalignment','center','verticalalignment','top')
xlabel('S')
ylabel('J dot')
%xlim([0 1.05*Tn2_0])
set(gca,'xtick',[0 Tn2_0])
set(gca,'xticklabels',{'0',''})
