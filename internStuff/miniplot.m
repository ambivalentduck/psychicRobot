function miniplot(x,y,xoff,marker,msize,rn,hh,hsize)

plot(xoff+rn,x-y,marker,'markersize',msize)

[~,h]=signrank(x,y);
if h
    plot(xoff,hh,'k*','markersize',hsize)
end



