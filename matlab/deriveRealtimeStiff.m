clc
clear all

%simplest example
H1=(M*s^2+B*s+K)/(M*s^2+B*s+K);



% syms t tf xi xf s
% syms M B K eps
% 
% ta=t/tf;
% xp=xi+(xf-xi)*(10*ta^3-15*ta^4+6*ta^5);
% %xp=xp*heaviside(-(t-tf))+xf*heaviside(t-tf)
% 
% Lxp=int(xp*exp(-s*t),t,0,tf)
% 
% simple(Lxp)
% pretty(Lxp)
% 
% ixp=inline(vectorize(xp))
% tv=0:.01:2;
% 
% x=ixp(tv,1.5,1,0);
% figure(1)
% clf
% hold on
% plot(tv,x)
% 
% 
% 
% H=(M*s^2+B*s+K)/(M*s^2+B*s+(K+eps));
% 
% FILT=H*Lxp;
% 
% filt=ilaplace(FILT)
% 
% ifilt=inline(vectorize(filt))
% 
% ifilt(1,1,1,s,tv,1.5,1,0)