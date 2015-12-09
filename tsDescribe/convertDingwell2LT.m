function convertDingwell2LT(X)

L=X(:,1);
T=X(:,2);

figure(1)
clf
subplot(2,1,1)
hold on
plotShiftedGam(T.^-2);
[h,pT]=jbtest(T)
[h,pTE]=jbtest(T.^-2)
pT-pTE
ylabel('Cumulative Frequecy')
xlabel('$\frac{1}{T^2}$','interpreter','latex')
%plotShiftedGam(T);
subplot(2,1,2)
hold on
plotShiftedGam((L./T).^2);
%plotShiftedGam(L);
[h,pL]=jbtest(L)
[h,pLE]=jbtest((L./T).^2)
pL-pLE
ylabel('Cumulative Frequecy')
xlabel('$\frac{L^2}{T^2}$','interpreter','latex')