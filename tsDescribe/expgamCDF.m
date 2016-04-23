function cdf=expgamCDF(x,nmin,nmean,theta)

lambda=nmean-nmin;
n=nmin+0:5*lambda;
npdf=exp(-n/lambda);
npdf=npdf/sum(npdf);

cdf=zeros(size(x));

for k=1:length(n)
    cdf=cdf+npdf(k)*gamcdf(x,n(k),theta);
end