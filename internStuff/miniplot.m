function [p,stats,mu]=miniplot(x,y,xoff,marker,msize,rn,hh,hsize)

plot(xoff+rn,x-y,marker,'markersize',msize)

%[p,h]=signrank(x,y);
[h,p,~,stats]=ttest(x,y);
if h
    plot(xoff,hh,'k*','markersize',hsize)
end

mu=mean(x-y);



