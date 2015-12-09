function out=getAlphaStatic(theta,omega,accel,l1,l2)
theta1=theta(1);
theta2=theta(2);

fJt1dx2=- l2.*cos(theta1 + theta2) - l1.*cos(theta1);
fJt1dy2=-l2.*cos(theta1 + theta2);
fJt2dx2=- l2.*sin(theta1 + theta2) - l1.*sin(theta1);
fJt2dy2=-l2.*sin(theta1 + theta2);
fJt1dxdy=-l2.*cos(theta1 + theta2);
fJt2dxdy=-l2.*sin(theta1 + theta2);
out=fJstatic(theta,l1,l2)\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);
