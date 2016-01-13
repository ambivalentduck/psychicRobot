clc
clear all

d=load('3dfreeexp1.dat');

t=d(:,1);
x=d(:,2:4);

figure(1)
clf
plot3(x(:,1),x(:,2),x(:,3))
xlabel('x')
ylabel('y')
zlabel('z')

v=zeros(size(x));
gT=gradient(t);
for k=1:3
    v(:,k)=gradient(x(:,k))./gT;
end

figure(2)
clf
subplot(2,1,1)
[f,h]=hist(x(:,3),15);
plot(h,log(f),'.')
xlabel('Height')
ylabel('ln Prob(height)')

z=x(:,3);
minz=min(z);
maxz=max(z);

nselect=15;
selectors=floor(15*(z-minz)/(maxz-minz))+1;
selectors(selectors==nselect+1)=nselect;
subplot(2,1,2)
plot(z,dot(v',v'),'.','markersize',.01)
xlabel('height')
ylabel('Speed^2')

figure(3)
clf
subplot(2,1,1)
vz=v(:,3);
speed=sqrt(dot(v',v'));
[fn,xn]=ecdf(speed(vz<0));
[fp,xp]=ecdf(speed(vz>0));
plot(xp,fp,'r',xn,fn,'b')
subplot(2,1,2)
[fn,xn]=ecdf(abs(vz(vz<0)));
[fp,xp]=ecdf(vz(vz>0));
plot(xp,fp,'r',xn,fn,'b')