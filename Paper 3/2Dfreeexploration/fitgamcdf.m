function cost=fitgamcdf(params)

global fitmecounts fitmebins

n=params(1);
T=params(2);

gcdf=gamcdf(fitmebins,n,T);
fit=gcdf(2:end)-gcdf(1:end-1);

cost=sum((fitmecounts-fit).^2);

