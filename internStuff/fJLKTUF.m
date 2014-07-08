function out=fJLKTUF(theta)
theta1=theta(1);
theta2=theta(2);
fJ11=- (381.*sin(theta1 + theta2))./1250 - (317.*sin(theta1))./1000;
fJ21=(381.*cos(theta1 + theta2))./1250 + (317.*cos(theta1))./1000;
fJ12=-(381.*sin(theta1 + theta2))./1250;
fJ22=(381.*cos(theta1 + theta2))./1250;
out=[fJ11,fJ12;fJ21,fJ22];
