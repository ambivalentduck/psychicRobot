function miniplot(x1,y1,x2,y2,msize,xl,rn)

hold on
plot(1+rn,x1-y1,'r.','markersize',msize)
plot(2+rn,x2-y2,'b.','markersize',msize)

yl=ylim;
a=.95;
asty=a*yl(2)+(1-a)*yl(1);

[~,h]=signrank(x1,y1);
if h
    plot(1,asty,'k*','markersize',msize)
end

[~,h]=signrank(x2,y2);
if h
    plot(2,asty,'k*','markersize',msize)
end

set(gca,'xtick',[1 2])
xlim(xl)
ylim(yl)


