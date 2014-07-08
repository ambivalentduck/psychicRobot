function out=getAlphaEYXMU(theta,omega,accel)
theta1=theta(1);
theta2=theta(2);
global fJ
fJt1dx2=- (9.*cos(theta1 + theta2))./25 - (7.*cos(theta1))./20;
fJt1dy2=-(9.*cos(theta1 + theta2))./25;
fJt2dx2=- (9.*sin(theta1 + theta2))./25 - (7.*sin(theta1))./20;
fJt2dy2=-(9.*sin(theta1 + theta2))./25;
fJt1dxdy=-(9.*cos(theta1 + theta2))./25;
fJt2dxdy=-(9.*sin(theta1 + theta2))./25;
out=fJ(theta)\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);
