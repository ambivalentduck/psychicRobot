function out=fJELWUY(theta)
theta1=theta(1);
theta2=theta(2);
fJ11=- (19.*sin(theta1 + theta2))./50 - (33.*sin(theta1))./100;
fJ21=(19.*cos(theta1 + theta2))./50 + (33.*cos(theta1))./100;
fJ12=-(19.*sin(theta1 + theta2))./50;
fJ22=(19.*cos(theta1 + theta2))./50;
out=[fJ11,fJ12;fJ21,fJ22];
