function out=rotateProgressError(x,x0,x1)

x=[x(:,1)-x0(1), x(:,2)-x0(2)];

xdiff=x1-x0;

theta=0-atan2(xdiff(2),xdiff(1));
ct=cos(theta);
st=sin(theta);
rotmat=[ct st;-st ct];

out=x*rotmat;
