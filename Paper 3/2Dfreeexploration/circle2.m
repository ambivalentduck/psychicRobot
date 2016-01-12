function h = circle2(x,y,r,alpha,color)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[1,1]);
set(h,'edgealpha',0,'alpha',alpha,'color',color);