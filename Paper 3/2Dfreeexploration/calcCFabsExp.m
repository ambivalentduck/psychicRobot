clc
clear all

syms t x real
syms A positive
syms N
assume(N>=2)
assumeAlso(N,'integer')
assumptions(N)

CFae=int(exp(i*t*x)*1/A*exp(-abs(x)/A),x,-inf,inf)
CFag=int(exp(i*t*x)*1/(gamma(N)*A^N)*abs(x)^(N-1)*exp(-abs(x)/A),x,-inf,inf)

% CF=2/(A^2*t^2 + 1)
%pdf=int(1/(2*pi)*exp(-i*t*x)*CF,t,0,inf)