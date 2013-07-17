clc
clear all
close all

load ../Data/output11.dat
load ../Data/input11.dat
load ../Data/11.mat

set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

val=sqrt(sum(output11(:,[5 6]).^2,2))+sqrt(sum(output11(:,[7 8]).^2,2));
f=find(val<.01);
figure(1)
subplot(2,2,1)
hold on
plot(output11(f,7),output11(f,9),'.')
plot(output11(f,8),output11(f,10),'r.')
ylabel('Force Sensor, N')
xlabel('Handle Acceleration, m/s^2')
legend('X','Y')
subplot(2,2,2)
filt=ones(1,1)/1;
hold on
plot(output11(f,1),output11(f,9),'b.',output11(f,1),output11(f,10),'r.',output11(f,1),filter(filt,1,output11(f,9)),'c-',output11(f,1),filter(filt,1,output11(f,10)),'m-','LineWidth',3)

kick=find(sum(abs(input11(max(1,output11(f,1)-1),[4 5])),2));
plot(output11(f(kick),1),output11(f(kick),10),'kx')
xlabel('Trial')
ylabel('Force Sensor Offset, N')
legend('X','Y','X Trend','Y Trend','Kick Came Just Before')

o=params.origin;
op=ikin(o);
e=[filter(filt,1,output11(f,9))';filter(filt,1,output11(f,10))'];
kp=SnMKpGainStruct(1.5);
kp=kp{1};
diff=f;
normf=f;
o=o';
for k=1:length(f)
    diff(k)=norm(o-fkin(kp\e(:,k)+op));
    normf(k)=norm(e(:,k));
end
subplot(2,2,3)
plot(output11(f,1),diff*100,'.')
xlabel('Trial')
ylabel('Estimated extraction offset due to force sensor drift, cm')

subplot(2,2,4)
plot(normf,diff*100,'.')
xlabel('Force Sensor Drift *Magnitude*, N')
ylabel('Extraction Error *Magnitude*, cm')