clc
clear all

syms L1 L2 C S positive
syms t real

td=t-0; %C1=0
ta=td/1+.5; %S1=1
K1=(30*ta^2-60*ta^3+30*ta^4)/1;

td=t-C;
ta=td/S+.5;
K2=(30*ta^2-60*ta^3+30*ta^4)/S;

assume(S>0)
assume(C>0)
assumeAlso((C-S/2)<.5)

extra=simplify(2*K1*K2)
rextra=simplify(int(extra,t,C-S/2,.5))
a=solve(diff(rextra,C^2),C^2)