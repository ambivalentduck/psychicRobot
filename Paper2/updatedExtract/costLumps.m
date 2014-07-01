function [cost,grad]=costLumps(in)

global nLumps

%Step 1 is to unpack in (indI, endF, xi, xf)
tspan=in(1:2*nLumps);
xi=reshape(in(2*nLumps+1:3*nLumps),nLumps,2);
xf=reshape(in(3*nLumps+1:4*nLumps),nLumps,2);

