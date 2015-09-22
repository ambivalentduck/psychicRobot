clc
clear all

t=linspace(0,.2,200);
dT=mean(diff(t));

figure(1)
clf
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

N=200;
nsubs=zeros(N,1);

for k=1:N
    
    [t1,t2,x1,x2]=makeSpaghetti;
    nsubs(k)=length(t1);
    Etotal=0*t;
    
    for k=1:length(t1)
        inds=find((t>=t1(k))&(t<=t2(k)));
        tc=(t1(k)+t2(k))/2;
        ts=t2(k)-t1(k);
        ta=((t(inds)-tc)+ts/2)/ts;
        Esub=(x2(k)-x1(k))*(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
        Etotal(inds)=Etotal(inds)+Esub.^2;
    end
    subplot(2,1,1)
    plot(t,sqrt(Etotal),'b-','linewidth',.001)
    subplot(2,1,2)
    plot(t,cumtrapz(sqrt(Etotal))*dT/2.8,'b-','linewidth',.001)
end
xlabel('Time, seconds')
ylabel('Progress Dimension Position, m')
ylim([0 .5])

subplot(2,1,1)
ylim([0 30])
ylabel('Speed, m/s')