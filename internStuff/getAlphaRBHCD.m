function out=getAlphaRBHCD(theta,omega,accel)
theta1=theta(1);
theta2=theta(2);
global fJ
fJt1dx2=- (381.*cos(theta1 + theta2))./1250 - (317.*cos(theta1))./1000;
fJt1dy2=-(381.*cos(theta1 + theta2))./1250;
fJt2dx2=- (381.*sin(theta1 + theta2))./1250 - (317.*sin(theta1))./1000;
fJt2dy2=-(381.*sin(theta1 + theta2))./1250;
fJt1dxdy=-(381.*cos(theta1 + theta2))./1250;
fJt2dxdy=-(381.*sin(theta1 + theta2))./1250;
out=fJ(theta)\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);
