clc
clear all

t=linspace(0,1,200);
dT=mean(diff(t));

figure(1)
clf
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

N=100;
nsubs=zeros(N,1);

for k=1:N
    
    [x1,x2,t1,t2]=makeSpaghetti;
    nsubs(k)=length(t1);
    Etotal=0*t;
    
    for k=1:length(t1)
        inds=find((t>=t1(k))&(t<=t2(k)));
        tc=(t1(k)+t2(k))/2;
        ts=t2(k)-t1(k);
        ta=((t(inds)-tc))/ts+.5;
        Esub=(x2(k)-x1(k))*(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
        Etotal(inds)=Etotal(inds)+Esub.^2;
    end
    speed=sqrt(Etotal);
    progress=cumtrapz(speed)*dT;
    subplot(2,1,1)
    plot(t,sqrt(Etotal),'b-','linewidth',.001)
    subplot(2,1,2)
    plot(progress,speed,'b-','linewidth',.001)
end
ylabel('Speed, m/s')
xlabel('Progress Dimension Position, m')
ylim([0 1])

subplot(2,1,1)
ylim([0 1])
ylabel('Speed, m/s')
xlabel('Time, s')
