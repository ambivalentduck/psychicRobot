function [cost,fit]=fitgamoffset(Tx)

global binedges bincounts nT


gcdf=gamcdf(binedges+Tx,nT(1),nT(2));
fit=gcdf(2:end)-gcdf(1:end-1);

cost=sum((bincounts-fit).^2);

