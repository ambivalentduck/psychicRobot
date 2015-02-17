function [cost,fit]=fitexpcdf(lambda)

global binedges bincounts

ecdf=expcdf(binedges,lambda);
fit=ecdf(2:end)-ecdf(1:end-1);
cost=sum((bincounts-fit).^2);