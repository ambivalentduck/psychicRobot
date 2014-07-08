function out=getAlphaTSKQF(theta,omega,accel)
theta1=theta(1);
theta2=theta(2);
global fJ
fJt1dx2=- (889.*cos(theta1 + theta2))./2500 - (127.*cos(theta1))./400;
fJt1dy2=-(889.*cos(theta1 + theta2))./2500;
fJt2dx2=- (889.*sin(theta1 + theta2))./2500 - (127.*sin(theta1))./400;
fJt2dy2=-(889.*sin(theta1 + theta2))./2500;
fJt1dxdy=-(889.*cos(theta1 + theta2))./2500;
fJt2dxdy=-(889.*sin(theta1 + theta2))./2500;
out=fJ(theta)\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);
