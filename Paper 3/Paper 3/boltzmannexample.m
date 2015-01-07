clc
clear all

figure(10)
clf

n=0:.01:1;

[Gx,Gy]=meshgrid(n,n);

Z=1-exp(-(Gx-.2).^2/.14-(Gy-.6).^2/.2)-exp(-(Gx-.8).^2/.08-(Gy-.3).^2/.2)

Fref=1;
Zref=0;

blah=0;

utemp=3;
ltemp=-2;
ntemp=100;

for T=[logspace(ltemp,utemp,ntemp/2) logspace(utemp,ltemp,ntemp/2)]
    C=Fref*exp((Zref-Z)/T);
    C=C/sum(C(:));
    surf(Gx,Gy,Z,C)
    blah=blah+1;
    title(['Temp=',num2str(T)])
    set(gca,'cameraposition',[1.5651   -1.6363    8.8248])
    F(blah)=getframe;
end
    