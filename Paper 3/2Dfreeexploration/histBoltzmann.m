function [R2,K,bins,counts,fitcounts]=histBoltzmann(X,lpctile,upctile,doplot)

upper=prctile(X,upctile);
lower=prctile(X,lpctile);

[counts,bins]=hist(X((X<=upper)&(X>=lower)),linspace(lower,upper,64));
counts=counts/sum(counts);
f=find(counts~=0);

logcounts=log(counts(f))';
nzbins=bins(f);

fitparams=[logcounts ones(size(logcounts))]\-nzbins';
fitcounts=exp((-bins-fitparams(2))/fitparams(1));
K=fitparams(1);

covfit=cov(logcounts,nzbins)/(std(logcounts)*std(nzbins));
R2=covfit(1,2)^2;

if nargin>2
    if doplot
        hold on
        bar(bins,counts)
        plot(bins,fitcounts,'r','linewidth',3)
    end
end

