function out=fJEYXMU(theta)
theta1=theta(1);
theta2=theta(2);
fJ11=- (9.*sin(theta1 + theta2))./25 - (7.*sin(theta1))./20;
fJ21=(9.*cos(theta1 + theta2))./25 + (7.*cos(theta1))./20;
fJ12=-(9.*sin(theta1 + theta2))./25;
fJ22=(9.*cos(theta1 + theta2))./25;
out=[fJ11,fJ12;fJ21,fJ22];
