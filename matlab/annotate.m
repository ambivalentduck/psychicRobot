function h=annotate(p,d,labs,colors,alength)

if nargin==1
    for k=1:length(p)
        delete(p(k).t)
        delete(p(k).a(1))
        delete(p(k).a(2))
    end
    h=[];
    return
end

latexscale=1;

if length(alength)==1
    alength=alength*ones(size(d,1));
end

for k=1:length(labs)
    x=p(k,1);
    y=p(k,2);
    l=-alength(k)*d(k,:)/norm(d(k,:));
    h(k).a=arrow([x y]+latexscale*l,[x y],colors(k,:),.3,2);
    
    if d(k,1)>0
        ha='right';
    elseif d(k,1)<0
        ha='left';
    else
        ha='center';
    end
    
    if d(k,2)>0
        va='top';
    elseif d(k,2)<0
        va='bottom';
    else
        va='middle';
    end
    
    h(k).t=text(x+l(1), y+l(2),labs{k},'horizontalalignment',ha,'Verticalalignment',va,'color',colors(k,:));
end