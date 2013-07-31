clc
clear all

out=load('../Data/output500.dat');
in=load('../Data/input_forcetest.dat');

reach1=find(out(:,1)==1);
onep=out(reach1,2:4);
onev=out(reach1,5:6);

figure(1)
clf
subplot(2,1,1)
plot3(onep(:,2),onep(:,3),onep(:,1))
still=find((vecmag([onep(:,2) - -.07349,onep(:,3)-.411])<.001)&(vecmag(onev)<.01));
subplot(2,1,2)
plot(gradient(still))

ostill=reach1(still);

figure(2)
clf
plot(out(ostill,9),out(ostill,10))
statsx=regstats(out(ostill,9),out(ostill,2),'linear',{'beta','rsquare'});
statsx.beta
statsx.rsquare
statsy=regstats(out(ostill,10),out(ostill,2),'linear',{'beta','rsquare'});
statsy.beta
statsy.rsquare

figure(3)
clf
hold on
aves=zeros(49,2);
contributing=zeros(49,1);
for k=2:49
    reach=find(out(:,1)==k);
    stillL=find((vecmag([out(reach,3)-.15,out(reach,4)-.5])<.005)&(vecmag(out(reach,5:6)<.01)));
    stillR=find((vecmag([out(reach,3)+.15,out(reach,4)-.5])<.005)&(vecmag(out(reach,5:6)<.01)));
    if length(stillL)<length(stillR)
        still=stillR;
    else
        still=stillL;
    end
    contributing(k-1)=length(still);
    force=out(reach(still),9:10);
    plot(force(:,1),force(:,2)+k*50,'.')
    aves(k-1,:)=mean(force);
end

figure(4)
clf
subplot(2,1,1)
plot(1:49,aves(:,1),'b',1:49,aves(:,2),'r')
legend('X','Y')
ylabel('Force, N')
subplot(2,1,2)
plot(1:49,abs(gradient(aves(:,1))),'b',1:49,abs(gradient(aves(:,2))),'r')
legend('X','Y')
xlabel('Trial')
ylabel('Abs(Gradient(Force)), N')

load ../Data/11.mat
set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)
o=params.origin;
op=ikin(o);
kp=SnMKpGainStruct(1.5);
kp=kp{1};
sserror=1:49;
o=o';
gaves=[gradient(aves(:,1)) gradient(aves(:,2))];
for k=1:49
    sserror(k)=norm(o-fkin(kp\gaves(k,:)'+op));
end

figure(5)
clf
plot(1:49,sserror*100)
ylabel('Steady-state error, cm')
xlabel('Trial')


figure(6)
clf
hold on
bar(1:49,log10(contributing))
xlabel('Trial')
ylabel('Log_{10}(number of contributing data points)')

