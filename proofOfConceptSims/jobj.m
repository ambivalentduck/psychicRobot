function ev=jobj(p)

global Tn2 deltaU L2

offset=p(1);
mu=p(2);

T=(mu*Tn2+offset).^-.5;
J=deltaU./T-L2.*T.^-3;
ev=-mean(J);
