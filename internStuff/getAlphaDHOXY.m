function out=getAlphaDHOXY(theta,omega,accel)
theta1=theta(1);
theta2=theta(2);
global fJ
fJt1dx2=- (19.*cos(theta1 + theta2))./50 - (9.*cos(theta1))./25;
fJt1dy2=-(19.*cos(theta1 + theta2))./50;
fJt2dx2=- (19.*sin(theta1 + theta2))./50 - (9.*sin(theta1))./25;
fJt2dy2=-(19.*sin(theta1 + theta2))./50;
fJt1dxdy=-(19.*cos(theta1 + theta2))./50;
fJt2dxdy=-(19.*sin(theta1 + theta2))./50;
out=fJ(theta)\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);
