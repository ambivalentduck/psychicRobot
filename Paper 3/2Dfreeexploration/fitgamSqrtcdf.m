function [cost,fit]=fitgamSqrtcdf(params)

global binedges bincounts

n=params(1);
T=params(2);

gcdf=gamcdf(sqrt(binedges),n,T);
fit=gcdf(2:end)-gcdf(1:end-1);

cost=sum((bincounts-fit).^2);

