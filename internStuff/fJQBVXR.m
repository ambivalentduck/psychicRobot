function out=fJQBVXR(theta)
theta1=theta(1);
theta2=theta(2);
fJ11=- (19.*sin(theta1 + theta2))./50 - (17.*sin(theta1))./50;
fJ21=(19.*cos(theta1 + theta2))./50 + (17.*cos(theta1))./50;
fJ12=-(19.*sin(theta1 + theta2))./50;
fJ22=(19.*cos(theta1 + theta2))./50;
out=[fJ11,fJ12;fJ21,fJ22];
