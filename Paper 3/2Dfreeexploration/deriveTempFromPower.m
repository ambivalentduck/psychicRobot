clc
clear all

syms P W x y omega real
syms W dt N positive

i=sqrt(-1);

pdf=(1/W)*exp(-P/W)

cdf=int(pdf,P,0,x)
%Answer is W*s, therefore s=1/W
invcdfboltz=solve(cdf-y,x)

mean=int(P*pdf,P,0,inf)

invfouriertrueboltz=int((1/W)*exp(-P/W)*exp(i*omega*P),P,0,inf)
recoveredtrueboltz=int(1/(2*pi)*exp(-i*omega*P)*invfouriertrueboltz,oemega)

invfourierboltz=int(.5*(1/W)*exp(-dt*abs(P)/W)*exp(i*omega*P),P,-inf,inf)

isprobdist=int(.5*(1/W)*exp(-dt*abs(P)/W),P,-inf,inf)

%Look at this over a 1 sec span
Nsolution=limit(invfourierboltz^(10/dt),dt,.3)

pdist=1/(2*pi)*int(exp(-i*omega*P)*Nsolution,omega,-inf,inf)

int(exp(-x)/x^N)