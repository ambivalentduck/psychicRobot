function ts=getTsMetric(x,t,x0,xf)

buffer=.1;

p=rotateProgressError(x,x0,xf);
onset=find(p(:,1)>buffer,1,'first');
offset=find(p(:,1)>(1-buffer),1,'first');
ts=t(offset)-t(onset);

function out=rotateProgressError(x,x0,x1)

x=[x(:,1)-x0(1), x(:,2)-x0(2)];

xdiff=x1-x0;

theta=0-atan2(xdiff(2),xdiff(1));
ct=cos(theta);
st=sin(theta);
rotmat=[ct st;-st ct];

out=x*rotmat;
out(:,1)=out(:,1)/norm(xdiff);