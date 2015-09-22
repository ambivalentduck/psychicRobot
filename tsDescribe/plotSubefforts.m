function plotSubefforts(t1,t2,x1,x2)

t=linspace(min(t1),max(t2),200);
figure(1)
clf
subplot(2,1,1)
hold on

Etotal=0*t;

for k=1:length(t1)
    color=.3+.4*rand(1,3);
    
    inds=find((t>=t1(k))&(t<=t2(k)));
    tc=(t1(k)+t2(k))/2;
    ts=t2(k)-t1(k);
    ta=((t(inds)-tc)+ts/2)/ts;
    Esub=(x2(k)-x1(k))*(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
    plot(t(inds),Esub,'-','color',color)
    
    Etotal(inds)=Etotal(inds)+Esub.^2;
end

plot(t,sqrt(Etotal),'k-')

subplot(2,1,2)
plot(t,cumtrapz(sqrt(Etotal)*mean(diff(t))))