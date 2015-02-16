function [cost,fit]=fitgamcdf(params)

global binedges bincounts

n=params(1);
T=params(2);

gcdf=gamcdf(binedges,n,T);
fit=gcdf(2:end)-gcdf(1:end-1);

cost=sum((bincounts-fit).^2);

