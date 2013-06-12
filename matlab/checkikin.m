clc
clear all

global l1 l2 x0
l1=.34;
l2=.33;
x0=[0 .8];

x=-.34:.01:.34;
y=0:.01:1;

[X,Y]=meshgrid(x,y);
Z=X;

for k=1:length(y)
    for kk=1:length(x)
        Z(k,kk)=norm([X(k,kk),Y(k,kk)]'-fkin(ikinAbs([X(k,kk),Y(k,kk)])));
    end
end

%Z=X.^2+Y;

figure(3)
clf
surf(X,Y,Z)